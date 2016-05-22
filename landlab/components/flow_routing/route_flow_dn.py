#! /usr/env/python

"""Calculate single-path (steepest direction) flow directions.

Given a ModelGrid, calculates single-path (steepest direction) flow directions,
drainage area, and (optionally) discharge.

The "dn" in the name means that this is a generalization of the D8 algorithm,
for a grid in which a node has N neighbors (N might happen to be 8, or not).
"""
# Created GT Nov 2013
# Modified to save data to grid directly, DEJH March 2014

from __future__ import print_function

import landlab
import warnings
from landlab.components.flow_routing import flow_direction_DN
from landlab.components.flow_accum import flow_accum_bw
from landlab import FieldError, Component
from landlab import ModelParameterDictionary
from landlab import RasterModelGrid, VoronoiDelaunayGrid  # for type tests
from landlab.utils.decorators import use_file_name_or_kwds
import numpy


class FlowRouter(Component):

    """Single-path (steepest direction) flow routing.

    This class implements single-path (steepest direction) flow routing, and
    calculates flow directions, drainage area, and (optionally) discharge.

    Note that because this router is generalizd across both regular and
    irregular grids, perimeter nodes can NEVER contribute to the accumulating
    flux, even if the gradients from them point inwards to the main body of
    the grid. This is because under Landlab definitions, perimeter nodes lack
    cells, so cannot accumulate any discharge.

    The primary method of this class is :func:`run_one_step`.

    Construction::

        FlowRouter(grid, method='D8', runoff_rate=None)

    Parameters
    ----------
    grid : ModelGrid
        A grid.
    method : {'D8', 'D4'}, optional
        Routing method ('D8' is the default). This keyword has no effect for a
        Voronoi-based grid.
    runoff_rate : float, optional (m/time)
        If provided, sets the (spatially constant) runoff rate. If a spatially
        variable runoff rate is desired, use the input field
        'water__unit_flux_in'. If both the field and argument are present at
        the time of initialization, runoff_rate will *overwrite* the field.
        If neither are set, defaults to spatially constant unit input.
    """

    _name = 'DNFlowRouter'

    _input_var_names = ('topographic__elevation',
                        'water__unit_flux_in',
                        )

    _output_var_names = ('drainage_area',
                         'flow__receiver_node',
                         'topographic__steepest_slope',
                         'water__discharge',
                         'flow__upstream_node_order',
                         'flow__link_to_receiver_node',
                         'flow__sink_flag',
                         )

    _var_units = {'topographic__elevation': 'm',
                  'water__unit_flux_in': 'm/s',
                  'drainage_area': 'm**2',
                  'flow__receiver_node': '-',
                  'topographic__steepest_slope': '-',
                  'water__discharge': 'm**3/s',
                  'flow__upstream_node_order': '-',
                  'flow__link_to_receiver_node': '-',
                  'flow__sink_flag': '-',
                  }

    _var_mapping = {'topographic__elevation': 'node',
                    'water__unit_flux_in': 'node',
                    'drainage_area': 'node',
                    'flow__receiver_node': 'node',
                    'topographic__steepest_slope': 'node',
                    'water__discharge': 'node',
                    'flow__upstream_node_order': 'node',
                    'flow__link_to_receiver_node': 'node',
                    'flow__sink_flag': 'node',
                    }

    _var_doc = {
        'topographic__elevation': 'Land surface topographic elevation',
        'water__unit_flux_in':
            ('External volume water per area per time input to each node ' +
             '(e.g., rainfall rate)'),
        'drainage_area':
            "Upstream accumulated surface area contributing to the node's "
            "discharge",
        'flow__receiver_node':
            'Node array of receivers (node that receives flow from current '
            'node)',
        'topographic__steepest_slope':
            'Node array of steepest *downhill* slopes',
        'water__discharge': 'Discharge of water through each node',
        'flow__upstream_node_order':
            'Node array containing downstream-to-upstream ordered list of '
            'node IDs',
        'flow__link_to_receiver_node':
            'ID of link downstream of each node, which carries the discharge',
        'flow__sink_flag': 'Boolean array, True at local lows',
    }

    @use_file_name_or_kwds
    def __init__(self, grid, method='D8', runoff_rate=None, **kwds):
        # We keep a local reference to the grid
        self._grid = grid
        if method in ('D8', 'D4', None):
            self.method = method
        else:
            raise AttributeError("method argument must be 'D8' or 'D4'")

        # set up the grid type testing
        self._is_raster = isinstance(self._grid, RasterModelGrid)
        if not self._is_raster:
            self.method = None

        # We'll also keep track of the active links; if raster, then these are
        # the "D8" links; otherwise, it's just activelinks
        if self._is_raster:
            dal, d8f, d8t = grid._d8_active_links()
            self._active_links = dal
            self._activelink_from = d8f
            self._activelink_to = d8t
            # needs modifying in the loop if D4 (now done)
        else:
            self._active_links = grid.active_links
            self._activelink_from = grid._activelink_fromnode
            self._activelink_to = grid._activelink_tonode

        # test input variables are present:
        grid.at_node['topographic__elevation']
        try:
            grid.at_node['water__unit_flux_in']
        except FieldError:
            grid.add_empty('node', 'water__unit_flux_in')
            if runoff_rate is None:
                grid.at_node['water__unit_flux_in'].fill(1.)
            else:
                grid.at_node['water__unit_flux_in'].fill(runoff_rate)
        else:
            if runoff_rate is not None:
                print ("FlowRouter found both the field " +
                       "'water__unit_flux_in' and a provided float or " +
                       "array for the runoff_rate argument. THE FIELD IS " +
                       "BEING OVERWRITTEN WITH THE SUPPLIED RUNOFF_RATE!")
                grid.at_node['water__unit_flux_in'].fill(runoff_rate)

        # perform a test (for politeness!) that the old name for the water_in
        # field is not present:
        try:
            grid.at_node['water__discharge_in']
        except FieldError:
            pass
        else:
            warnings.warn("This component formerly took 'water__discharge" +
                          "_in' as an input field. However, this field is " +
                          "now named 'water__unit_flux_in'. You are still " +
                          "using a field with the old name. Please update " +
                          "your code if you intended the FlowRouter to use " +
                          "that field.", DeprecationWarning)

        # Keep track of the following variables:
        #   - drainage area at each node
        #   - receiver of each node
        self.drainage_area = grid.add_zeros('drainage_area', at='node',
                                            noclobber=False)
        self.receiver = grid.add_zeros('flow__receiver_node', at='node',
                                       noclobber=False, dtype=int)
        self.steepest_slope = grid.add_zeros('topographic__steepest_' +
                                             'slope', at='node',
                                             noclobber=False)
        self.discharges = grid.add_zeros('water__discharge', at='node',
                                         noclobber=False)
        self.upstream_ordered_nodes = grid.add_zeros('flow__upstream_node_order',
                                                     at='node', dtype=int,
                                                     noclobber=False)
        self.links_to_receiver = grid.add_zeros('flow__link_to_receiver_node',
                                                at='node', dtype=int,
                                                noclobber=False)
        grid.add_zeros('flow__sink_flag', at='node', dtype=int, noclobber=False)

    def route_flow(self, **kwds):
        """Route surface-water flow over a landscape.

        Routes surface-water flow by (1) assigning to each node a single
        drainage direction, and then (2) adding up the number of nodes that
        contribute flow to each node on the grid (including the node itself).

        Stores as ModelGrid fields:

        -  Node array of receivers (nodes that receive flow), or ITS OWN ID if
           there is no receiver: *'flow__receiver_node'*
        -  Node array of drainage areas: *'drainage_area'*
        -  Node array of discharges: *'water__discharge'*
        -  Node array of steepest downhill slopes:
           *'topographic__steepest_slope'*
        -  Node array containing downstream-to-upstream ordered list of node
           IDs: *'flow__upstream_node_order'*
        -  Node array containing ID of link that leads from each node to its
           receiver, or BAD_INDEX_VALUE if no link:
           *'flow__link_to_receiver_node'*
        -  Boolean node array of all local lows: *'flow__sink_flag'*

        Returns
        -------
        ModelGrid
            The modified grid object

        Examples
        --------
        >>> import numpy as np
        >>> from landlab import RasterModelGrid
        >>> from landlab.components.flow_routing import FlowRouter
        >>> mg = RasterModelGrid((5, 4), spacing=(1, 1))
        >>> elev = np.array([0.,  0.,  0., 0.,
        ...                  0., 21., 10., 0.,
        ...                  0., 31., 20., 0.,
        ...                  0., 32., 30., 0.,
        ...                  0.,  0.,  0., 0.])
        >>> _ = mg.add_field('node','topographic__elevation', elev)
        >>> mg.set_closed_boundaries_at_grid_edges(True, True, True, False)
        >>> fr = FlowRouter(mg)
        >>> mg = fr.route_flow()
        >>> mg.at_node['flow__receiver_node'] # doctest: +NORMALIZE_WHITESPACE
        array([  0,  1,  2,  3,
                 4,  1,  2,  7,
                 8,  6,  6, 11,
                12, 10, 10, 15,
                16, 17, 18, 19])
        >>> mg.at_node['drainage_area'] # doctest: +NORMALIZE_WHITESPACE
        array([ 0.,  1.,  5.,  0.,
                0.,  1.,  5.,  0.,
                0.,  1.,  3.,  0.,
                0.,  1.,  1.,  0.,
                0.,  0.,  0.,  0.])

        Now let's change the cell area (100.) and the runoff rates:

        >>> mg = RasterModelGrid((5, 4), spacing=(10., 10))

        Put the data back into the new grid.

        >>> _ = mg.add_field('node','topographic__elevation', elev)
        >>> mg.set_closed_boundaries_at_grid_edges(True, True, True, False)
        >>> fr = FlowRouter(mg)
        >>> runoff_rate = np.arange(mg.number_of_nodes)
        >>> _ = mg.add_field('node', 'water__unit_flux_in', runoff_rate,
        ...                  noclobber=False)
        >>> mg = fr.route_flow()
        >>> mg.at_node['water__discharge'] # doctest: +NORMALIZE_WHITESPACE
        array([    0.,   500.,  5200.,     0.,
                   0.,   500.,  5200.,     0.,
                   0.,   900.,  3700.,     0.,
                   0.,  1300.,  1400.,     0.,
                   0.,     0.,     0.,     0.])

        """
        # this retained for back compatibility - method now set in __init__.
        if 'method' in kwds:
            warnings.warn("'method' should be set at initialization now. " +
                          "Please update your code.", DeprecationWarning)
            # raise NameError
            if kwds['method'] not in ('D8', 'D4'):
                raise ValueError('method not understood ({method})'.format(
                    method=method))
            else:
                self.method = kwds['method']
            if not self._is_raster:
                self.method = None

        # if elevs is not provided, default to stored grid values, which must
        # be provided as grid
        elevs = self._grid['node']['topographic__elevation']

        node_cell_area = self._grid.cell_area_at_node.copy()
        node_cell_area[self._grid.closed_boundary_nodes] = 0.
        # closed cells can't contribute

        # Calculate the downhill-positive slopes at the d8 active links
        if self.method == 'D8':
            link_slope = - self._grid._calculate_gradients_at_d8_active_links(
                elevs)
        else:
            link_slope = - self._grid.calc_grad_of_active_link(
                elevs)

        # Find the baselevel nodes
        (baselevel_nodes, ) = numpy.where(
            numpy.logical_or(self._grid.status_at_node == 1,
                             self._grid.status_at_node == 2))

        # Calculate flow directions
        if self.method == 'D4':
            num_d4_active = self._grid.number_of_active_links  # only d4
            receiver, steepest_slope, sink, recvr_link = \
                flow_direction_DN.flow_directions(elevs, self._active_links,
                                         self._activelink_from[:num_d4_active],
                                         self._activelink_to[:num_d4_active],
                                         link_slope,
                                         grid=self._grid,
                                         baselevel_nodes=baselevel_nodes)
        else:  # Voronoi or D8
            receiver, steepest_slope, sink, recvr_link = \
                flow_direction_DN.flow_directions(elevs, self._active_links,
                                     self._activelink_from,
                                     self._activelink_to, link_slope,
                                     grid=self._grid,
                                     baselevel_nodes=baselevel_nodes)

        # TODO: either need a way to calculate and return the *length* of the
        # flow links, OR the caller has to handle the raster / non-raster case.

        # Calculate drainage area, discharge, and ...
        a, q, s = flow_accum_bw.flow_accumulation(
            receiver, sink, node_cell_area=node_cell_area,
            runoff_rate=self._grid.at_node['water__unit_flux_in'])

        # added DEJH March 2014:
        # store the generated data in the grid
        self._grid['node']['drainage_area'][:] = a
        self._grid['node']['flow__receiver_node'][:] = receiver
        self._grid['node']['topographic__steepest_slope'][:] = steepest_slope
        self._grid['node']['water__discharge'][:] = q
        self._grid['node']['flow__upstream_node_order'][:] = s
        self._grid['node']['flow__link_to_receiver_node'][:] = recvr_link
        self._grid['node']['flow__sink_flag'][:] = numpy.zeros_like(receiver,
                                                               dtype=bool)
        self._grid['node']['flow__sink_flag'][sink] = True

        return self._grid

    def run_one_step(self, **kwds):
        """Route surface-water flow over a landscape.

        Routes surface-water flow by (1) assigning to each node a single
        drainage direction, and then (2) adding up the number of nodes that
        contribute flow to each node on the grid (including the node itself).

        This is the fully standardized run method for this component. It
        differs from :func:`route_flow` in that it has a standardized name,
        and does not return anything.

        Examples
        --------
        >>> import numpy as np
        >>> from landlab import RasterModelGrid
        >>> from landlab.components.flow_routing import FlowRouter
        >>> mg = RasterModelGrid((5, 4), spacing=(1, 1))
        >>> elev = np.array([0.,  0.,  0., 0.,
        ...                  0., 21., 10., 0.,
        ...                  0., 31., 20., 0.,
        ...                  0., 32., 30., 0.,
        ...                  0.,  0.,  0., 0.])
        >>> _ = mg.add_field('node','topographic__elevation', elev)
        >>> mg.set_closed_boundaries_at_grid_edges(True, True, True, False)
        >>> fr = FlowRouter(mg)
        >>> fr.run_one_step()
        >>> mg.at_node['flow__receiver_node'] # doctest: +NORMALIZE_WHITESPACE
        array([  0,  1,  2,  3,
                 4,  1,  2,  7,
                 8,  6,  6, 11,
                12, 10, 10, 15,
                16, 17, 18, 19])
        >>> mg.at_node['drainage_area'] # doctest: +NORMALIZE_WHITESPACE
        array([ 0.,  1.,  5.,  0.,
                0.,  1.,  5.,  0.,
                0.,  1.,  3.,  0.,
                0.,  1.,  1.,  0.,
                0.,  0.,  0.,  0.])

        Now let's change the cell area (100.) and the runoff rates:

        >>> mg = RasterModelGrid((5, 4), spacing=(10., 10))

        Put the data back into the new grid.

        >>> _ = mg.add_field('node','topographic__elevation', elev)
        >>> mg.set_closed_boundaries_at_grid_edges(True, True, True, False)
        >>> fr = FlowRouter(mg)
        >>> runoff_rate = np.arange(mg.number_of_nodes)
        >>> _ = mg.add_field('node', 'water__unit_flux_in', runoff_rate,
        ...                  noclobber=False)
        >>> fr.run_one_step()
        >>> mg.at_node['water__discharge'] # doctest: +NORMALIZE_WHITESPACE
        array([    0.,   500.,  5200.,     0.,
                   0.,   500.,  5200.,     0.,
                   0.,   900.,  3700.,     0.,
                   0.,  1300.,  1400.,     0.,
                   0.,     0.,     0.,     0.])
        """
        self.route_flow(**kwds)

    @property
    def node_drainage_area(self):
        return self._grid['node']['drainage_area']

    @property
    def node_receiving_flow(self):
        return self._grid['node']['flow__receiver_node']

    @property
    def node_steepest_slope(self):
        return self._grid['node']['topographic__steepest_slope']

    @property
    def node_water_discharge(self):
        return self._grid['node']['water__discharge']

    @property
    def node_order_upstream(self):
        return self._grid['node']['flow__upstream_node_order']

    @property
    def link_to_flow_receiving_node(self):
        return self._grid['node']['flow__link_to_receiver_node']


if __name__ == '__main__':
    import doctest
    doctest.testmod()

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
from landlab.components.flow_routing import flow_direction_DN
from landlab.components.flow_accum import flow_accum_bw
from landlab import FieldError, Component
from landlab import ModelParameterDictionary
from landlab import RasterModelGrid, VoronoiDelaunayGrid  # for type tests
import numpy


#output_suppression_flag = True

class FlowRouter(Component):

    """Single-path (steepest direction) flow routing.

    This class implements single-path (steepest direction) flow routing, and
    calculates flow directions, drainage area, and (optionally) discharge.

    It initializes with a reference to a ModelGrid of any kind. Optionally, it
    can also take *input_params*, the string which is the name of a text input
    file. The input file can optionally contain 'runoff_rate', a float giving
    a (spatially constant) runoff rate. This is equivalent to the optional
    input field 'water__volume_flux_in', and will override it if both are set.
    If neither are set, value will default to 1.

    Note that because this router is generalizd across both regular and
    irregular grids, perimeter nodes can NEVER contribute to the accumulating
    flux, even if the gradients from them point inwards to the main body of
    the grid. This is because under Landlab definitions, perimeter nodes lack
    cells, so cannot accumulate any discharge.

    The primary method of this class is :func:`route_flow`.
    """

    _name = 'DNFlowRouter'

    _input_var_names = set(['topographic__elevation',
                            'water__volume_flux_in',
                            ])

    _output_var_names = set(['drainage_area',
                             'flow_receiver',
                             'topographic__steepest_slope',
                             'water__volume_flux',
                             'upstream_node_order',
                             'links_to_flow_receiver',
                             'flow_sinks',
                             ])

    _var_units = {'topographic__elevation' : 'm',
                  'water__volume_flux_in' : 'm**3/s',
                  'drainage_area' : 'm**2',
                  'flow_receiver' : '-',
                  'topographic__steepest_slope' : '-',
                  'water__volume_flux' : 'm**3/s',
                  'upstream_node_order' : '-',
                  'links_to_flow_receiver' : '-',
                  'flow_sinks' : '-',
                  }

    _var_mapping = {'topographic__elevation' : 'node',
                    'water__volume_flux_in' : 'node',
                    'drainage_area' : 'node',
                    'flow_receiver' : 'node',
                    'topographic__steepest_slope' : 'node',
                    'water__volume_flux' : 'node',
                    'upstream_node_order' : 'node',
                    'links_to_flow_receiver' : 'node',
                    'flow_sinks' : 'node',
                    }

    _var_doc = {
        'topographic__elevation': 'Land surface topographic elevation',
        'water__volume_flux_in':
            'External volume water input to each node (e.g., rainfall)',
        'drainage_area':
            "Upstream accumulated surface area contributing to the node's "
            "discharge",
        'flow_receiver':
            'Node array of receivers (node that receives flow from current '
            'node)',
        'topographic__steepest_slope':
            'Node array of steepest *downhill* slopes',
        'water__volume_flux': 'Discharge of water through each node',
        'upstream_node_order':
            'Node array containing downstream-to-upstream ordered list of '
            'node IDs',
        'links_to_flow_receiver':
            'ID of link downstream of each node, which carries the discharge',
        'flow_sinks': 'Boolean array, True at local lows',
    }

    def __init__(self, model_grid, input_params=None):
        # We keep a local reference to the grid
        self._grid = model_grid
        self.value_field = 'topographic__elevation'

        # set up the grid type testing
        self._is_raster = isinstance(self._grid, RasterModelGrid)

        if input_params:
            if type(input_params) == str:
                input_dict = ModelParameterDictionary(input_params)
            else:
                assert type(input_params) == dict
                input_dict = input_params

        # We'll also keep track of the active links; if raster, then these are
        # the "D8" links; otherwise, it's just activelinks
        if self._is_raster:
            dal, d8f, d8t = model_grid.d8_active_links()
            self._active_links = dal
            self._activelink_from = d8f
            self._activelink_to = d8t
            # needs modifying in the loop if D4 (now done)
        else:
            self._active_links = model_grid.active_links
            self._activelink_from = model_grid.activelink_fromnode
            self._activelink_to = model_grid.activelink_tonode

        # test input variables are present:
        model_grid.at_node['topographic__elevation']
        try:
            model_grid.at_node['water__volume_flux_in']
            self.field_for_runoff = True
        except FieldError:
            self.field_for_runoff = False
            # build the input array into the grid. This is important in case
            # variable values appear during model run
            model_grid.add_ones('node', 'water__volume_flux_in')

        if not self.field_for_runoff:
            try:
                model_grid.at_node['water__volume_flux_in'].fill(
                    input_dict['runoff_rate'])
            except (KeyError, UnboundLocalError):
                model_grid.at_node['water__volume_flux_in'].fill(1.)
        else:
            try:
                model_grid.at_node['water__volume_flux_in'].fill(
                    input_dict['runoff_rate'])
            except (KeyError, UnboundLocalError):
                pass
            else:
                print("WARNING: Both a field and input parameter are "
                      "available for runoff value. Was this intentional?? "
                      "Taking the input parameter value...")

        # Keep track of the following variables:
        #   - drainage area at each node
        #   - receiver of each node
        self.drainage_area = model_grid.add_zeros('node', 'drainage_area')
        self.receiver = model_grid.create_node_array_zeros('flow_receiver')
        self.steepest_slope = model_grid.create_node_array_zeros(
            'topographic__steepest_slope')
        self.discharges = model_grid.create_node_array_zeros(
            'water__volume_flux')
        self.upstream_ordered_nodes = model_grid.create_node_array_zeros(
            'upstream_node_order')
        self.links_to_receiver = model_grid.create_node_array_zeros(
            'links_to_flow_receiver')

    def route_flow(self, method='D8'):
        """Route surface-water flow over a landscape.

        Routes surface-water flow by (1) assigning to each node a single
        drainage direction, and then (2) adding up the number of nodes that
        contribute flow to each node on the grid (including the node itself).

        Stores as ModelGrid fields:
        - Node array of receivers (nodes that receive flow):
              *'flow_receiver'*
        - Node array of drainage areas: *'drainage_area'*
        - Node array of discharges: *'water__volume_flux'*
        - Node array of steepest downhill slopes:
          *'topographic__steepest_slope'*
        - Node array containing downstream-to-upstream ordered list of node
          IDs: *'upstream_node_order'*
        - Node array containing ID of link that leads from each node to its
          receiver (or ITS OWN ID if there is no receiver):
          *'links_to_flow_receiver'*
        - Boolean node array of all local lows: *'flow_sinks'*

        Parameters
        ----------
        method : {'D8', 'D4'}, optional
            Routing method ('D8' is the default). This keyword has no effect
            for a Voronoi-based grid.

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
        >>> mg.set_closed_boundaries_at_grid_edges(False, True, True, True)
        >>> fr = FlowRouter(mg)
        >>> mg = fr.route_flow()
        >>> mg.at_node['flow_receiver'] # doctest: +NORMALIZE_WHITESPACE
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
        >>> mg.set_closed_boundaries_at_grid_edges(False, True, True, True)
        >>> fr = FlowRouter(mg)
        >>> runoff_rate = np.arange(mg.number_of_nodes)
        >>> _ = mg.add_field('node', 'water__volume_flux_in', runoff_rate)
        >>> mg = fr.route_flow()
        >>> mg.at_node['water__volume_flux'] # doctest: +NORMALIZE_WHITESPACE
        array([    0.,   500.,  5200.,     0.,
                   0.,   500.,  5200.,     0.,
                   0.,   900.,  3700.,     0.,
                   0.,  1300.,  1400.,     0.,
                   0.,     0.,     0.,     0.])

        """
        if method not in ('D8', 'D4'):
            raise ValueError('method not understood ({method})'.format(
                method=method))
        if not self._is_raster:
            method = None

        # if elevs is not provided, default to stored grid values, which must
        # be provided as grid
        elevs = self._grid['node'][self.value_field]

        node_cell_area = self._grid.cell_area_at_node.copy()
        node_cell_area[self._grid.closed_boundary_nodes] = 0.
        # closed cells can't contribute


        # Calculate the downhill-positive slopes at the d8 active links
        if method == 'D8':
            link_slope = - self._grid.calculate_gradients_at_d8_active_links(
                elevs)
        else:
            link_slope = - self._grid.calculate_gradients_at_active_links(
                elevs)

        # Find the baselevel nodes
        (baselevel_nodes, ) = numpy.where(
            numpy.logical_or(self._grid.status_at_node == 1,
                             self._grid.status_at_node == 2))

        # Calculate flow directions
        if method == 'D4':
            num_d4_active = self._grid.number_of_active_links  # only d4
            receiver, steepest_slope, sink, recvr_link  = \
                flow_direction_DN.flow_directions(elevs, self._active_links,
                                         self._activelink_from[:num_d4_active],
                                         self._activelink_to[:num_d4_active],
                                         link_slope,
                                         grid=self._grid,
                                         baselevel_nodes=baselevel_nodes)
        else:  # Voronoi or D8
            receiver, steepest_slope, sink, recvr_link  = \
                flow_direction_DN.flow_directions(elevs, self._active_links,
                                     self._activelink_from,
                                     self._activelink_to, link_slope,
                                     grid=self._grid,
                                     baselevel_nodes=baselevel_nodes)
        #############grid=None???

        # TODO: either need a way to calculate and return the *length* of the
        # flow links, OR the caller has to handle the raster / non-raster case.

        #print 'sinks:', sink

        # Calculate drainage area, discharge, and ...
        a, q, s = flow_accum_bw.flow_accumulation(
            receiver, sink, node_cell_area=node_cell_area,
            runoff_rate=self._grid.at_node['water__volume_flux_in'])

        #added DEJH March 2014:
        #store the generated data in the grid
        self._grid['node']['drainage_area'] = a
        self._grid['node']['flow_receiver'] = receiver
        self._grid['node']['topographic__steepest_slope'] = steepest_slope
        self._grid['node']['water__volume_flux'] = q
        self._grid['node']['upstream_node_order'] = s
        self._grid['node']['links_to_flow_receiver'] = recvr_link
        self._grid['node']['flow_sinks'] = numpy.zeros_like(receiver,
                                                            dtype=bool)
        self._grid['node']['flow_sinks'][sink] = True

        return self._grid

    @property
    def node_drainage_area(self):
        return self._grid['node']['drainage_area']

    @property
    def node_receiving_flow(self):
        return self._grid['node']['flow_receiver']

    @property
    def node_steepest_slope(self):
        return self._grid['node']['topographic__steepest_slope']

    @property
    def node_water_discharge(self):
        return self._grid['node']['water__volume_flux']

    @property
    def node_order_upstream(self):
        return self._grid['node']['upstream_node_order']

    @property
    def link_to_flow_receiving_node(self):
        return self._grid['node']['links_to_flow_receiver']


if __name__ == '__main__':
    import doctest
    doctest.testmod()

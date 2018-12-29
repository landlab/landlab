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

#
# import landlab
import warnings

# from landlab.components.flow_director import flow_direction_DN
# from landlab.components.flow_accum import flow_accum_bw
# from landlab import FIXED_VALUE_BOUNDARY, FIXED_GRADIENT_BOUNDARY
from landlab import VoronoiDelaunayGrid  # for type tests
from landlab.components.flow_accum.flow_accumulator import FlowAccumulator
from landlab.utils.decorators import use_file_name_or_kwds


class FlowRouter(FlowAccumulator):

    """Single-path (steepest direction) flow routing.

    **AS OF v1.5.3 THIS COMPONENT HAS BEEN DEPRECATED. USE THE FlowAccumulator
    INSTEAD. THIS COMPONENT WILL BE REMOVED IN v2.0**

    This class implements single-path (steepest direction) flow routing, and
    calculates flow directions, drainage area, and (optionally) discharge.

    Note that because this router is generalizd across both regular and
    irregular grids, perimeter nodes can NEVER contribute to the accumulating
    flux, even if the gradients from them point inwards to the main body of
    the grid. This is because under Landlab definitions, perimeter nodes lack
    cells, so cannot accumulate any discharge.

    The primary method of this class is :func:`run_one_step`.
    """

    _name = "DNFlowRouter"

    _input_var_names = ("topographic__elevation", "water__unit_flux_in")

    _output_var_names = (
        "drainage_area",
        "flow__receiver_node",
        "topographic__steepest_slope",
        "surface_water__discharge",
        "flow__upstream_node_order",
        "flow__link_to_receiver_node",
        "flow__sink_flag",
    )

    _var_units = {
        "topographic__elevation": "m",
        "water__unit_flux_in": "m/s",
        "drainage_area": "m**2",
        "flow__receiver_node": "-",
        "topographic__steepest_slope": "-",
        "surface_water__discharge": "m**3/s",
        "flow__upstream_node_order": "-",
        "flow__link_to_receiver_node": "-",
        "flow__sink_flag": "-",
    }

    _var_mapping = {
        "topographic__elevation": "node",
        "water__unit_flux_in": "node",
        "drainage_area": "node",
        "flow__receiver_node": "node",
        "topographic__steepest_slope": "node",
        "surface_water__discharge": "node",
        "flow__upstream_node_order": "node",
        "flow__link_to_receiver_node": "node",
        "flow__sink_flag": "node",
    }

    _var_doc = {
        "topographic__elevation": "Land surface topographic elevation",
        "water__unit_flux_in": (
            "External volume water per area per time input to each node "
            + "(e.g., rainfall rate)"
        ),
        "drainage_area": "Upstream accumulated surface area contributing to the node's "
        "discharge",
        "flow__receiver_node": "Node array of receivers (node that receives flow from current "
        "node)",
        "topographic__steepest_slope": "Node array of steepest *downhill* slopes",
        "surface_water__discharge": "Discharge of water through each node",
        "flow__upstream_node_order": "Node array containing downstream-to-upstream ordered list of "
        "node IDs",
        "flow__link_to_receiver_node": "ID of link downstream of each node, which carries the discharge",
        "flow__sink_flag": "Boolean array, True at local lows",
    }

    @use_file_name_or_kwds
    def __init__(self, grid, method="D8", runoff_rate=None, **kwds):
        """Initialize FlowDirector.

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
        msg = (
            "FlowRouter has been deprecated as of Landlab v1.5.2 and will be "
            "removed in v2.0. Use FlowAccumulator instead."
        )
        warnings.warn(msg, DeprecationWarning)
        self._is_Voroni = isinstance(grid, VoronoiDelaunayGrid)
        self._grid = grid
        if "method" in kwds:
            warnings.warn(
                "'method' should be set at initialization now. "
                + "Please update your code.",
                DeprecationWarning,
            )
            # raise NameError
            if kwds["method"] not in ("D8", "D4"):
                raise ValueError(
                    "method not understood ({method})".format(method=method)
                )

        if method == "D4" or self._is_Voroni:
            flow_director = "Steepest"
        else:
            flow_director = "D8"

        super(FlowRouter, self).__init__(
            grid,
            surface="topographic__elevation",
            flow_director=flow_director,
            runoff_rate=runoff_rate,
        )

    def _test_for_method_change(self, **kwds):
        """Provides backwards compatability for method keyword.

        The original flow router allowed the method to be specified as a
        keyword argument to run_one_step. Under the new flow accumulator
        framework this requires resetting which flow director is used.
        """

        # this retained for back compatibility - method now set in __init__.
        if "method" in kwds:

            method = kwds.pop("method")
            warnings.warn(
                "'method' should be set at initialization now. "
                + "Please update your code.",
                DeprecationWarning,
            )
            # raise NameError
            if method not in ("D8", "D4"):
                raise ValueError(
                    "method not understood ({method})".format(method=method)
                )
            else:
                self.method = method
            if not self._is_raster:
                self.method = None

            if method == "D4" or self._is_Voroni:
                flow_director = "Steepest"
            else:
                flow_director = "D8"

            self._add_director(flow_director)

    def route_flow(self, **kwds):
        """Route surface-water flow over a landscape.

        Routes surface-water flow by (1) assigning to each node a single
        drainage direction, and then (2) adding up the number of nodes that
        contribute flow to each node on the grid (including the node itself).

        Stores as ModelGrid fields:

        -  Node array of receivers (nodes that receive flow), or ITS OWN ID if
           there is no receiver: *'flow__receiver_node'*
        -  Node array of drainage areas: *'drainage_area'*
        -  Node array of discharges: *'surface_water__discharge'*
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
        >>> import pytest
        >>> import numpy as np
        >>> from landlab import RasterModelGrid
        >>> from landlab.components.flow_routing import FlowRouter
        >>> mg = RasterModelGrid((5, 4), xy_spacing=(1, 1))
        >>> elev = np.array([0.,  0.,  0., 0.,
        ...                  0., 21., 10., 0.,
        ...                  0., 31., 20., 0.,
        ...                  0., 32., 30., 0.,
        ...                  0.,  0.,  0., 0.])
        >>> _ = mg.add_field('node','topographic__elevation', elev)
        >>> mg.set_closed_boundaries_at_grid_edges(True, True, True, False)
        >>> with pytest.deprecated_call():
        ...    fr = FlowRouter(mg)
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

        >>> mg = RasterModelGrid((5, 4), xy_spacing=(10., 10))

        Put the data back into the new grid.

        >>> _ = mg.add_field('node','topographic__elevation', elev)
        >>> mg.set_closed_boundaries_at_grid_edges(True, True, True, False)
        >>> with pytest.deprecated_call():
        ...     fr = FlowRouter(mg)
        >>> runoff_rate = np.arange(mg.number_of_nodes)
        >>> _ = mg.add_field('node', 'water__unit_flux_in', runoff_rate,
        ...                  noclobber=False)
        >>> mg = fr.route_flow()
        >>> mg.at_node['surface_water__discharge'] # doctest: +NORMALIZE_WHITESPACE
        array([    0.,   500.,  5200.,     0.,
                   0.,   500.,  5200.,     0.,
                   0.,   900.,  3700.,     0.,
                   0.,  1300.,  1400.,     0.,
                   0.,     0.,     0.,     0.])
        """
        self._test_for_method_change(**kwds)
        self.accumulate_flow()

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
        >>> import pytest
        >>> import numpy as np
        >>> from landlab import RasterModelGrid
        >>> from landlab.components.flow_routing import FlowRouter
        >>> mg = RasterModelGrid((5, 4), xy_spacing=(1, 1))
        >>> elev = np.array([0.,  0.,  0., 0.,
        ...                  0., 21., 10., 0.,
        ...                  0., 31., 20., 0.,
        ...                  0., 32., 30., 0.,
        ...                  0.,  0.,  0., 0.])
        >>> _ = mg.add_field('node','topographic__elevation', elev)
        >>> mg.set_closed_boundaries_at_grid_edges(True, True, True, False)
        >>> with pytest.deprecated_call():
        ...    fr = FlowRouter(mg)
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

        The default behavior of FlowRouter is to use the D8 method. Next we
        will examine the alternative case of the D4 method that does not
        consider diagonal links bewtween nodes.

        >>> mg = RasterModelGrid((5, 4), xy_spacing=(1, 1))
        >>> elev = np.array([0.,  0.,  0., 0.,
        ...                  0., 21., 10., 0.,
        ...                  0., 31., 20., 0.,
        ...                  0., 32., 30., 0.,
        ...                  0.,  0.,  0., 0.])
        >>> _ = mg.add_field('node','topographic__elevation', elev)
        >>> mg.set_closed_boundaries_at_grid_edges(True, True, True, False)
        >>> with pytest.deprecated_call():
        ...     fr = FlowRouter(mg, method = 'D4')
        >>> fr.run_one_step()
        >>> mg.at_node['flow__receiver_node'] # doctest: +NORMALIZE_WHITESPACE
        array([ 0,  1,  2,  3,
                4,  1,  2,  7,
                8, 10,  6, 11,
               12, 14, 10, 15,
               16, 17, 18, 19])
        >>> mg.at_node['drainage_area'] # doctest: +NORMALIZE_WHITESPACE
        array([ 0.,  1.,  5.,  0.,
                0.,  1.,  5.,  0.,
                0.,  1.,  4.,  0.,
                0.,  1.,  2.,  0.,
                0.,  0.,  0.,  0.])

        The flow router can also work on irregular grids.  For the example we
        will use a Hexagonal Model Grid, a special type of Voroni Grid that has
        regularly spaced hexagonal cells. We will also set the dx spacing such
        that each cell has an area of one.

        >>> from landlab import HexModelGrid
        >>> dx=(2./(3.**0.5))**0.5
        >>> mg = HexModelGrid(5,3, dx)
        >>> _ = mg.add_field('topographic__elevation', mg.node_x + np.round(mg.node_y), at = 'node')
        >>> with pytest.deprecated_call():
        ...    fr = FlowRouter(mg)
        >>> fr.run_one_step()
        >>> mg.at_node['flow__receiver_node'] # doctest: +NORMALIZE_WHITESPACE
        array([ 0,  1,  2,
                3,  0,  1,  6,
                7,  3,  4,  5, 11,
               12,  8,  9, 15,
               16, 17, 18])
        >>> mg.at_node['drainage_area'] # doctest: +NORMALIZE_WHITESPACE
        array([ 3.,  2.,  0.,
                2.,  3.,  2.,  0.,
                0.,  2.,  2.,  1.,  0.,
                0., 1.,  1.,  0.,
                0.,  0.,  0.])

        Now let's return to the first example and change the cell area (100.)
        and the runoff rates:

        >>> mg = RasterModelGrid((5, 4), xy_spacing=(10., 10))

        Put the data back into the new grid.

        >>> _ = mg.add_field('node','topographic__elevation', elev)
        >>> mg.set_closed_boundaries_at_grid_edges(True, True, True, False)
        >>> with pytest.deprecated_call():
        ...    fr = FlowRouter(mg)
        >>> runoff_rate = np.arange(mg.number_of_nodes)
        >>> _ = mg.add_field('node', 'water__unit_flux_in', runoff_rate,
        ...                  noclobber=False)
        >>> fr.run_one_step()
        >>> mg.at_node['surface_water__discharge'] # doctest: +NORMALIZE_WHITESPACE
        array([    0.,   500.,  5200.,     0.,
                   0.,   500.,  5200.,     0.,
                   0.,   900.,  3700.,     0.,
                   0.,  1300.,  1400.,     0.,
                   0.,     0.,     0.,     0.])
        """

        self._test_for_method_change(**kwds)
        self.accumulate_flow()

    @property
    def node_drainage_area(self):
        return self._grid["node"]["drainage_area"]

    @property
    def node_receiving_flow(self):
        return self._grid["node"]["flow__receiver_node"]

    @property
    def node_steepest_slope(self):
        return self._grid["node"]["topographic__steepest_slope"]

    @property
    def node_water_discharge(self):
        return self._grid["node"]["surface_water__discharge"]

    @property
    def node_order_upstream(self):
        return self._grid["node"]["flow__upstream_node_order"]

    @property
    def link_to_flow_receiving_node(self):
        return self._grid["node"]["flow__link_to_receiver_node"]


if __name__ == "__main__":
    import doctest

    doctest.testmod()

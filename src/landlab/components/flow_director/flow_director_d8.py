#! /usr/env/python

"""
flow_director_d8.py: provides the component FlowDirectorsD8.

This components finds the steepest single-path steepest descent flow
directions and considers diagonal links between nodes on a raster grid. It is
not implemented for irregular grids. For a method that works for irregular
grids and does not consider diagonal links for rasters, use
FlowDirectorSteepest instead.
"""

import numpy

from landlab import LinkStatus
from landlab.components.flow_director import flow_direction_DN
from landlab.components.flow_director.flow_director_to_one import _FlowDirectorToOne


class FlowDirectorD8(_FlowDirectorToOne):
    """Single-path (steepest direction) flow direction with diagonals on
    rasters.

    Single-path (steepest direction) flow direction finding on raster grids
    by the D8 method. This method considers flow on all eight links such that
    flow is possible on orthogonal and on diagonal links.

    The method that considers only orthogonal links (D4 method) for raster
    grids is FlowDirectorSteepest.

    This method is not implemented for Voroni grids, use
    FlowDirectorSteepest instead.

    Stores as ModelGrid fields:

    -  Node array of receivers (nodes that receive flow), or ITS OWN ID if
       there is no receiver: *'flow__receiver_node'*
    -  Node array of steepest downhill slopes:
       *'topographic__steepest_slope'*
    -  Node array containing ID of link that leads from each node to its
       receiver, or BAD_INDEX_VALUE if no link:
       *'flow__link_to_receiver_node'*
    -  Boolean node array of all local lows: *'flow__sink_flag'*

    The primary method of this class is :func:`run_one_step`.

    Examples
    --------
    >>> import numpy as np
    >>> from landlab import RasterModelGrid
    >>> from landlab.components import FlowDirectorD8
    >>> mg = RasterModelGrid((3, 3), xy_spacing=(1, 1))
    >>> mg.set_closed_boundaries_at_grid_edges(True, True, True, False)
    >>> _ = mg.add_field(
    ...     "topographic__elevation",
    ...     mg.node_x + mg.node_y,
    ...     at="node",
    ... )
    >>> fd = FlowDirectorD8(mg, "topographic__elevation")
    >>> fd.surface_values
    array([0., 1., 2., 1., 2., 3., 2., 3., 4.])
    >>> fd.run_one_step()
    >>> mg.at_node["flow__receiver_node"]
    array([0, 1, 2, 3, 0, 5, 6, 7, 8])
    >>> mg.at_node["topographic__steepest_slope"]
    array([0.        , 0.        , 0.        , 0.        , 1.41421356,
           0.        , 0.        , 0.        , 0.        ])
    >>> mg.at_node["flow__link_to_receiver_node"]
    array([-1, -1, -1, -1, 12, -1, -1, -1, -1])
    >>> mg.at_node["flow__sink_flag"].astype(int)
    array([1, 1, 1, 1, 0, 1, 1, 1, 1])
    >>> mg_2 = RasterModelGrid((5, 4), xy_spacing=(1, 1))
    >>> topographic__elevation = [
    ...     [0.0, 0.0, 0.0, 0.0],
    ...     [0.0, 21.0, 10.0, 0.0],
    ...     [0.0, 31.0, 20.0, 0.0],
    ...     [0.0, 32.0, 30.0, 0.0],
    ...     [0.0, 0.0, 0.0, 0.0],
    ... ]
    >>> _ = mg_2.add_field(
    ...     "topographic__elevation",
    ...     topographic__elevation,
    ...     at="node",
    ... )
    >>> mg_2.set_closed_boundaries_at_grid_edges(True, True, True, False)
    >>> fd_2 = FlowDirectorD8(mg_2)
    >>> fd_2.run_one_step()
    >>> mg_2.at_node["flow__receiver_node"].reshape(mg_2.shape)
    array([[ 0,  1,  2,  3],
           [ 4,  1,  2,  7],
           [ 8,  6,  6, 11],
           [12, 10, 10, 15],
           [16, 17, 18, 19]])

    The flow directors also have the ability to return the flow receiver nodes

    >>> receiver = fd.direct_flow()
    >>> receiver
    array([0, 1, 2,
           3, 0, 5,
           6, 7, 8])

    References
    ----------
    **Required Software Citation(s) Specific to this Component**

    None Listed

    **Additional References**

    O'Callaghan, J., Mark, D. (1984). The extraction of drainage networks from
    digital elevation data. Computer Vision, Graphics, and Image Processing
    28(3), 323 - 344. https://dx.doi.org/10.1016/s0734-189x(84)80011-0

    """

    _name = "FlowDirectorD8"

    _info = {
        "flow__link_to_receiver_node": {
            "dtype": int,
            "intent": "out",
            "optional": False,
            "units": "-",
            "mapping": "node",
            "doc": "ID of link downstream of each node, which carries the discharge",
        },
        "flow__receiver_node": {
            "dtype": int,
            "intent": "out",
            "optional": False,
            "units": "-",
            "mapping": "node",
            "doc": "Node array of receivers (node that receives flow from current node)",
        },
        "flow__sink_flag": {
            "dtype": bool,
            "intent": "out",
            "optional": False,
            "units": "-",
            "mapping": "node",
            "doc": "Boolean array, True at local lows",
        },
        "topographic__elevation": {
            "dtype": float,
            "intent": "in",
            "optional": True,
            "units": "m",
            "mapping": "node",
            "doc": "Land surface topographic elevation",
        },
        "topographic__steepest_slope": {
            "dtype": float,
            "intent": "out",
            "optional": False,
            "units": "-",
            "mapping": "node",
            "doc": "The steepest *downhill* slope",
        },
    }

    def __init__(self, grid, surface="topographic__elevation"):
        """
        Parameters
        ----------
        grid : ModelGrid
            A grid of type RasterModelGrid.
        surface : field name at node or array of length node, optional
            The surface to direct flow across, default is field at node:
            topographic__elevation,.
        """
        self._method = "D8"
        super().__init__(grid, surface)
        try:
            self._grid.nodes_at_d8
        except AttributeError:
            self._is_Voroni = True
        else:
            self._is_Voroni = False
        if self._is_Voroni:
            raise NotImplementedError(
                "FlowDirectorD8 not implemented for"
                "irregular grids, use"
                "FlowDirectorSteepest"
            )

        self.updated_boundary_conditions()

    def updated_boundary_conditions(self):
        """Method to update FlowDirectorD8 when boundary conditions change.

        Call this if boundary conditions on the grid are updated after
        the component is instantiated.
        """
        self._active_links = numpy.arange(self._grid.number_of_d8)
        nodes_at_d8 = self._grid.nodes_at_d8[self._active_links]
        self._activelink_tail = nodes_at_d8[:, 0]
        self._activelink_head = nodes_at_d8[:, 1]

    def run_one_step(self):
        """Find flow directions and save to the model grid.

        run_one_step() checks for updated boundary conditions, calculates
        slopes on links, finds baselevel nodes based on the status at node,
        calculates flow directions, and saves results to the grid.

        an alternative to direct_flow() is direct_flow() which does the same
        things but also returns the receiver nodes not return values.
        """
        self.direct_flow()

    def direct_flow(self):
        """Find flow directions, save to the model grid, and return receivers.

        direct_flow() checks for updated boundary conditions, calculates
        slopes on links, finds baselevel nodes based on the status at node,
        calculates flow directions, saves results to the grid, and returns a
        at-node array  of receiver nodes. This array is stored in the grid at:
        grid['node']['flow__receiver_node']

        an alternative to direct_flow() is run_one_step() which does the same
        things but also returns a at-node array  of receiver nodes. This array
        is stored in the grid at:
        grid['node']['flow__receiver_node']
        """
        self._check_updated_bc()

        # update the surface, if it was provided as a model grid field.
        self._changed_surface()

        # step 1. Calculate link slopes.
        link_slope = -self._grid.calc_grad_at_d8(self._surface_values)
        link_slope[self._grid.status_at_d8 != LinkStatus.ACTIVE] = 0

        # Step 2. Find and save base level nodes.
        (baselevel_nodes,) = numpy.where(
            numpy.logical_or(
                self._grid.status_at_node == self._grid.BC_NODE_IS_FIXED_VALUE,
                self._grid.status_at_node == self._grid.BC_NODE_IS_FIXED_GRADIENT,
            )
        )

        # Calculate flow directions by D8 method
        receiver, steepest_slope, sink, recvr_link = flow_direction_DN.flow_directions(
            self._surface_values,
            self._active_links,
            self._activelink_tail,
            self._activelink_head,
            link_slope,
            grid=self._grid,
            baselevel_nodes=baselevel_nodes,
        )
        # Save the four ouputs of this component.
        self._grid["node"]["flow__receiver_node"][:] = receiver
        self._grid["node"]["topographic__steepest_slope"][:] = steepest_slope
        self._grid["node"]["flow__link_to_receiver_node"][:] = recvr_link
        self._grid["node"]["flow__sink_flag"][:] = numpy.zeros_like(
            receiver, dtype=bool
        )
        self._grid["node"]["flow__sink_flag"][sink] = True

        return receiver

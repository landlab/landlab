"""Find depressions on a topographic surface.

.. codeauthor:: gtucker, DEJH (Flow routing)
"""

# Routing by DEJH, Oct 15.


import numpy as np

from ...core.model_component import Component
from ...core.utils import as_id_array
from ...field import FieldError
from ...grid import RasterModelGrid
from ..flow_accum import flow_accum_bw
from .cfuncs import find_lowest_node_on_lake_perimeter_c
from .floodstatus import FloodStatus

_UNFLOODED = FloodStatus._UNFLOODED
_CURRENT_LAKE = FloodStatus._CURRENT_LAKE
_FLOODED = FloodStatus._FLOODED
_PIT = FloodStatus._PIT

use_cfuncs = True


class DepressionFinderAndRouter(Component):
    """Find depressions on a topographic surface.

    This component identifies depressions in a topographic surface, finds an
    outlet for each depression.  If directed to do so (default True), and the
    component is able to find existing routing fields output from the
    'route_flow_dn' component, it will then modify the drainage directions and
    accumulations already stored in the grid to route flow across these
    depressions.

    Note that in general properties of this class named "depression" identify
    each individual pit in the topography, including those that will merge
    once the fill is performed. Those named "lake" return the unique lakes
    created by the fill, and are probably the properties most users will
    want.

    Note also that the structure of drainage within the lakes is not
    guaranteed, and in particular, may not be symmetrical even if your
    boundary conditions are.
    However, the outputs from the lake will all still be correct.

    Note the routing part of this component may not yet be fully compatible with
    irregular grids.

    The prinary method of this class is
    *map_depressions()*.

    Examples
    --------
    Route flow across a depression in a sloped surface.

    >>> from landlab import RasterModelGrid
    >>> from landlab.components import FlowAccumulator, DepressionFinderAndRouter
    >>> mg = RasterModelGrid((7, 7), xy_spacing=0.5)
    >>> z = mg.add_field("topographic__elevation", mg.node_x.copy(), at="node")
    >>> z += 0.01 * mg.node_y
    >>> mg.at_node["topographic__elevation"].reshape(mg.shape)[2:5, 2:5] *= 0.1
    >>> fr = FlowAccumulator(mg, flow_director="D8")
    >>> fr.run_one_step()  # the flow "gets stuck" in the hole
    >>> mg.at_node["flow__receiver_node"].reshape(mg.shape)
    array([[ 0,  1,  2,  3,  4,  5,  6],
           [ 7,  7, 16, 17, 18, 18, 13],
           [14, 14, 16, 16, 17, 18, 20],
           [21, 21, 16, 23, 24, 25, 27],
           [28, 28, 23, 30, 31, 32, 34],
           [35, 35, 30, 31, 32, 32, 41],
           [42, 43, 44, 45, 46, 47, 48]])
    >>> mg.at_node["drainage_area"].reshape(mg.shape)
    array([[0.  , 0.  , 0.  , 0.  , 0.  , 0.  , 0.  ],
           [0.25, 0.25, 0.25, 0.25, 0.25, 0.25, 0.  ],
           [0.25, 0.25, 5.  , 1.5 , 1.  , 0.25, 0.  ],
           [0.25, 0.25, 3.  , 0.75, 0.5 , 0.25, 0.  ],
           [0.25, 0.25, 2.  , 1.5 , 1.  , 0.25, 0.  ],
           [0.25, 0.25, 0.25, 0.25, 0.25, 0.25, 0.  ],
           [0.  , 0.  , 0.  , 0.  , 0.  , 0.  , 0.  ]])
    >>> df = DepressionFinderAndRouter(mg)
    >>> df.map_depressions()  # reroute_flow defaults to True
    >>> mg.at_node["flow__receiver_node"].reshape(mg.shape)
    array([[ 0,  1,  2,  3,  4,  5,  6],
           [ 7,  7, 16, 17, 18, 18, 13],
           [14, 14,  8, 16, 17, 18, 20],
           [21, 21, 16, 16, 24, 25, 27],
           [28, 28, 23, 24, 24, 32, 34],
           [35, 35, 30, 31, 32, 32, 41],
           [42, 43, 44, 45, 46, 47, 48]])
    >>> mg.at_node["drainage_area"].reshape(mg.shape)
    array([[0.  , 0.  , 0.  , 0.  , 0.  , 0.  , 0.  ],
           [5.25, 5.25, 0.25, 0.25, 0.25, 0.25, 0.  ],
           [0.25, 0.25, 5.  , 1.5 , 1.  , 0.25, 0.  ],
           [0.25, 0.25, 0.75, 2.25, 0.5 , 0.25, 0.  ],
           [0.25, 0.25, 0.5 , 0.5 , 1.  , 0.25, 0.  ],
           [0.25, 0.25, 0.25, 0.25, 0.25, 0.25, 0.  ],
           [0.  , 0.  , 0.  , 0.  , 0.  , 0.  , 0.  ]])
    >>> df.lake_at_node.reshape(mg.shape)
    array([[False, False, False, False, False, False, False],
           [False, False, False, False, False, False, False],
           [False, False,  True,  True,  True, False, False],
           [False, False,  True,  True,  True, False, False],
           [False, False,  True,  True,  True, False, False],
           [False, False, False, False, False, False, False],
           [False, False, False, False, False, False, False]])
    >>> df.lake_map.reshape(mg.shape)
    array([[-1, -1, -1, -1, -1, -1, -1],
           [-1, -1, -1, -1, -1, -1, -1],
           [-1, -1, 16, 16, 16, -1, -1],
           [-1, -1, 16, 16, 16, -1, -1],
           [-1, -1, 16, 16, 16, -1, -1],
           [-1, -1, -1, -1, -1, -1, -1],
           [-1, -1, -1, -1, -1, -1, -1]])
    >>> df.lake_codes  # a unique code for each lake present on the grid
    array([16])
    >>> df.lake_outlets  # the outlet node of each lake in lake_codes
    array([8])
    >>> df.lake_areas  # the area of each lake in lake_codes
    array([2.25])

    Because ``rereoute_flow`` defaults to ``True``, the flow connectivity fields
    created by the :py:class:`~landlab.components.flow_accum.FlowAccumulator`
    will have now been modified to route flow over the depressions in the
    surface. The topogrphy itself is not modified.

    References
    ----------
    **Required Software Citation(s) Specific to this Component**

    None Listed

    **Additional References**

    Tucker, G. E., Lancaster, S. T., Gasparini, N. M., and Bras, R. L.: The
    Channel-Hillslope Integrated Landscape Development Model (CHILD), in:
    Landscape Erosion and Evolution Modeling, Springer US, Boston, MA, USA,
    349–388, 2001.

    """

    _name = "DepressionFinderAndRouter"

    _unit_agnostic = True

    _info = {
        "depression__depth": {
            "dtype": float,
            "intent": "out",
            "optional": False,
            "units": "m",
            "mapping": "node",
            "doc": "Depth of depression below its spillway point",
        },
        "depression__outlet_node": {
            "dtype": int,
            "intent": "out",
            "optional": False,
            "units": "-",
            "mapping": "node",
            "doc": (
                "If a depression, the id of the outlet node for that depression, "
                "otherwise grid.BAD_INDEX"
            ),
        },
        "flood_status_code": {
            "dtype": int,
            "intent": "out",
            "optional": False,
            "units": "-",
            "mapping": "node",
            "doc": "Map of flood status (_PIT, _CURRENT_LAKE, _UNFLOODED, or _FLOODED).",
        },
        "is_pit": {
            "dtype": bool,
            "intent": "out",
            "optional": False,
            "units": "-",
            "mapping": "node",
            "doc": "Boolean flag indicating whether a node is a pit.",
        },
        "topographic__elevation": {
            "dtype": float,
            "intent": "in",
            "optional": False,
            "units": "m",
            "mapping": "node",
            "doc": "Land surface topographic elevation",
        },
    }

    def __init__(self, grid, routing="D8", pits="flow__sink_flag", reroute_flow=True):
        """Create a DepressionFinderAndRouter.

        Constructor assigns a copy of the grid, sets the current time, and
        calls the initialize method.

        Parameters
        ----------
        grid : RasterModelGrid
            A landlab RasterModelGrid.
        routing : str
            If grid is a raster type, controls whether lake connectivity can
            occur on diagonals ('D8', default), or only orthogonally ('D4').
            Has no effect if grid is not a raster.
        pits : array or str or None, optional
            If a field name, the boolean field containing True where pits.
            If an array, either a boolean array of nodes of the pits, or an
            array of pit node IDs. It does not matter whether or not open
            boundary nodes are flagged as pits; they are never treated as such.
            Default is 'flow__sink_flag', the pit field output from the
            :py:mod:`FlowDirectors <landlab.components.flow_director>`.
        reroute_flow : bool, optional
            If True (default), and the component detects the output fields in
            the grid produced by the FlowAccumulator component, this component
            will modify the existing flow fields to route the flow across the
            lake surface(s).
        """
        super().__init__(grid)

        self._bc_set_code = self._grid.bc_set_code

        self._user_supplied_pits = pits
        self._reroute_flow = reroute_flow

        if routing != "D8":
            assert routing == "D4"
        self._routing = routing

        if isinstance(grid, RasterModelGrid) and (routing == "D8"):
            self._D8 = True
            self._num_nbrs = 8
            self._diag_link_length = np.sqrt(grid.dx**2 + grid.dy**2)
        else:
            self._D8 = False  # useful shorthand for thia test we do a lot
            if isinstance(grid, RasterModelGrid):
                self._num_nbrs = 4
            else:
                self._num_nbrs = self._grid.links_at_node.shape[1]

        if "flow__receiver_node" in self._grid.at_node and self._grid.at_node[
            "flow__receiver_node"
        ].size != self._grid.size("node"):
            raise NotImplementedError(
                "A route-to-multiple flow director has been "
                "run on this grid. The depression finder is "
                "not compatible with the grid anymore. Use "
                "DepressionFinderAndRouter with reroute_flow=True "
                "only with route-to-one methods. If using this "
                "component with such a flow directing method is desired "
                "please open a GitHub Issue/"
            )

        # Make sure the grid includes elevation data.
        self._elev = self._grid.at_node["topographic__elevation"]

        # Create output variables.
        #
        # Note that we initialize depression
        # outlet ID to self._grid.BAD_INDEX (which is a major clue!)
        self._depression_depth = self._grid.add_zeros(
            "depression__depth", at="node", clobber=True
        )
        self._depression_outlet_map = self._grid.add_zeros(
            "depression__outlet_node", at="node", dtype=int, clobber=True
        )
        self._depression_outlet_map += self._grid.BAD_INDEX

        # Later on, we'll need a number that's guaranteed to be larger than the
        # highest elevation in the grid.
        self._BIG_ELEV = 1.0e99

        self.updated_boundary_conditions()

        self._lake_outlets = []  # a list of each unique lake outlet
        # ^note this is nlakes-long

        self._is_pit = self._grid.add_ones(
            "is_pit", at="node", dtype=bool, clobber=True
        )
        self._flood_status = self._grid.add_zeros(
            "flood_status_code", at="node", dtype=int, clobber=True
        )
        self._lake_map = np.empty(self._grid.number_of_nodes, dtype=int)
        self._lake_map.fill(self._grid.BAD_INDEX)

    def updated_boundary_conditions(self):
        """Call this if boundary conditions on the grid are updated after the
        component is instantiated."""
        try:
            dx = self._grid.dx
            dy = self._grid.dy
        except AttributeError:
            pass
        # We'll also need a handy copy of the node neighbor lists
        # TODO: presently, this grid method seems to only exist for Raster
        # grids. We need it for *all* grids!
        self._node_nbrs = self._grid.active_adjacent_nodes_at_node
        if self._D8:
            diag_nbrs = self._grid.diagonal_adjacent_nodes_at_node.copy()
            # remove the inactive nodes:
            diag_nbrs[
                self._grid.status_at_node[diag_nbrs] == self._grid.BC_NODE_IS_CLOSED
            ] = -1
            self._node_nbrs = np.concatenate((self._node_nbrs, diag_nbrs), 1)
            self._link_lengths = np.empty(8, dtype=float)
            self._link_lengths[0] = dx
            self._link_lengths[2] = dx
            self._link_lengths[1] = dy
            self._link_lengths[3] = dy
            self._link_lengths[4:].fill(np.sqrt(dx * dx + dy * dy))
        elif isinstance(self._grid, RasterModelGrid) and (self._routing == "D4"):
            self._link_lengths = np.empty(4, dtype=float)
            self._link_lengths[0] = dx
            self._link_lengths[2] = dx
            self._link_lengths[1] = dy
            self._link_lengths[3] = dy
        else:
            self._link_lengths = self._grid.length_of_link

    @property
    def is_pit(self):
        """At node array indicating whether the node is a pit or not."""
        return self._is_pit

    @property
    def number_of_pits(self):
        """The number of pits on the grid."""
        return self._number_of_pits

    @property
    def pit_node_ids(self):
        """Node IDs of grid nodes identified as pits."""
        return self._pit_node_ids

    @property
    def flood_status(self):
        """Map of flood status (_PIT, _CURRENT_LAKE, _UNFLOODED, or
        _FLOODED)."""
        return self._flood_status

    @property
    def receivers(self):
        """At node array indicating which node receives flow."""
        return self._receivers

    @receivers.setter
    def receivers(self, receivers):
        self._receivers = receivers

    @property
    def depression_depth(self):
        """At node array of depression depths."""
        return self._depression_depth

    @property
    def depression_outlet_map(self):
        """At node array indicating the node-id of the depression outlet."""
        return self._depression_outlet_map

    def _find_pits(self):
        """Locate local depressions ("pits") in a gridded elevation field.

        Notes
        -----
        **Uses**:

        * ``self._elev``
        * ``self._grid``

        **Creates**:

        * ``is_pit`` (node array of booleans): Flag indicating whether
          the node is a pit.
        * ``number_of_pits`` (int): Number of pits found.
        * ``pit_node_ids`` (node array of ints): IDs of the nodes that
          are pits

        A node is defined as being a pit if and only if:

        1. All neighboring core nodes have equal or greater elevation, and
        2. Any neighboring open boundary nodes have a greater elevation.

        The algorithm starts off assuming that all core nodes are pits. We then
        loop through all active links. For each link, if one node is higher
        than the other, the higher one cannot be a pit, so we flag it False.
        We also look at cases in which an active link's nodes have equal
        elevations. If one is an open boundary, then the other must be a core
        node, and we declare the latter not to be a pit (via rule 2 above).
        """
        # Create the is_pit array, with all core nodes initialized to True and
        # all boundary nodes initialized to False.
        self._is_pit.fill(True)
        self._is_pit[self._grid.boundary_nodes] = False

        # Loop over all active links: if one of a link's two nodes is higher
        # than the other, the higher one is not a pit. Also, if they have
        # equal elevations and one is an open boundary, the other is not a pit.
        act_links = self._grid.active_links
        h_orth = self._grid.node_at_link_head[act_links]
        t_orth = self._grid.node_at_link_tail[act_links]

        # These two lines assign the False flag to any node that is higher
        # than its partner on the other end of its link
        self._is_pit[h_orth[np.where(self._elev[h_orth] > self._elev[t_orth])[0]]] = (
            False
        )
        self._is_pit[t_orth[np.where(self._elev[t_orth] > self._elev[h_orth])[0]]] = (
            False
        )

        # If we have a raster grid, handle the diagonal active links too
        # (At the moment, their data structure is a bit different)
        # TODO: update the diagonal link data structures
        # DEJH doesn't understand why this can't be vectorized as above...
        if self._D8:
            for h, t in self._grid.nodes_at_diagonal[self._grid.active_diagonals]:
                if self._elev[h] > self._elev[t]:
                    self._is_pit[h] = False
                elif self._elev[t] > self._elev[h]:
                    self._is_pit[t] = False
                elif self._elev[h] == self._elev[t]:
                    if (
                        self._grid.status_at_node[h]
                        == self._grid.BC_NODE_IS_FIXED_VALUE
                    ):
                        self._is_pit[t] = False
                    elif (
                        self._grid.status_at_node[t]
                        == self._grid.BC_NODE_IS_FIXED_VALUE
                    ):
                        self._is_pit[h] = False

        # Record the number of pits and the IDs of pit nodes.
        self._number_of_pits = np.count_nonzero(self._is_pit)
        self._pit_node_ids = as_id_array(np.where(self._is_pit)[0])

    def _links_and_nbrs_at_node(self, the_node):
        """Compile and return arrays with IDs of neighbor links and nodes.

        If D8 Raster, returns *diag_nbrs* containing the diagonal neighbors;
        otherwise, *diag_nbrs = None*.

        Examples
        --------
        >>> from landlab import RasterModelGrid
        >>> from landlab.components import DepressionFinderAndRouter
        >>> rg = RasterModelGrid((3, 3))
        >>> z = rg.add_zeros("topographic__elevation", at="node")
        >>> z[4] = 2.0
        >>> df = DepressionFinderAndRouter(rg, routing="D4")
        >>> (links, nbrs, dnbrs) = df._links_and_nbrs_at_node(4)
        >>> links
        array([6, 8, 5, 3])
        >>> nbrs
        array([5, 7, 3, 1])
        >>> dnbrs
        >>> df = DepressionFinderAndRouter(rg, routing="D8")
        >>> (links, nbrs, dnbrs) = df._links_and_nbrs_at_node(4)
        >>> links
        array([6, 8, 5, 3])
        >>> nbrs
        array([5, 7, 3, 1])
        >>> dnbrs
        array([8, 6, 0, 2])
        """

        # Get the neighboring links (and, if applicable, the diagonals)
        links = self._grid.links_at_node[the_node]
        nbrs = self._grid.adjacent_nodes_at_node[the_node]
        if self._D8:
            diag_nbrs = self._grid.diagonal_adjacent_nodes_at_node[the_node]
        else:
            diag_nbrs = None

        return links, nbrs, diag_nbrs

    def assign_outlet_receiver(self, outlet_node):
        """Find drainage direction for outlet_node that does not flow into its
        own lake.

        Examples
        --------
        >>> import numpy as np
        >>> from landlab.components import DepressionFinderAndRouter
        >>> from landlab.components import FlowAccumulator
        >>> from landlab import RasterModelGrid
        >>> rg = RasterModelGrid((7, 7))
        >>> rg.status_at_node[rg.nodes_at_right_edge] = rg.BC_NODE_IS_CLOSED
        >>> z = rg.add_zeros("topographic__elevation", at="node")
        >>> z[:] = rg.x_of_node + 0.01 * rg.y_of_node
        >>> lake_nodes = np.array([10, 16, 17, 18, 24, 32, 33, 38, 40])
        >>> z[lake_nodes] *= 0.1
        >>> fr = FlowAccumulator(rg, flow_director="D4")
        >>> fr.run_one_step()
        >>> rg.at_node["flow__receiver_node"].reshape(rg.shape)
        array([[ 0,  1,  2,  3,  4,  5,  6],
               [ 7,  7, 16, 10, 10, 11, 13],
               [14, 14, 16, 16, 17, 18, 20],
               [21, 21, 16, 17, 24, 33, 27],
               [28, 28, 29, 24, 32, 32, 34],
               [35, 35, 38, 38, 38, 33, 41],
               [42, 43, 44, 45, 46, 47, 48]])
        >>> df = DepressionFinderAndRouter(rg, routing="D4")
        >>> df.map_depressions()
        >>> rg.at_node["flow__receiver_node"].reshape(rg.shape)
        array([[ 0,  1,  2,  3,  4,  5,  6],
               [ 7,  7, 16, 17, 10, 11, 13],
               [14, 14, 15, 16, 17, 18, 20],
               [21, 21, 16, 17, 24, 33, 27],
               [28, 28, 29, 38, 31, 32, 34],
               [35, 35, 36, 37, 38, 33, 41],
               [42, 43, 44, 45, 46, 47, 48]])
        >>> fr = FlowAccumulator(rg, flow_director="D8")
        >>> fr.run_one_step()
        >>> rg.at_node["flow__receiver_node"].reshape(rg.shape)
        array([[ 0,  1,  2,  3,  4,  5,  6],
               [ 7,  7, 16, 16, 10, 18, 13],
               [14, 14, 16, 16, 17, 18, 20],
               [21, 21, 16, 16, 24, 33, 27],
               [28, 28, 24, 24, 24, 32, 34],
               [35, 35, 38, 38, 38, 32, 41],
               [42, 43, 44, 45, 46, 47, 48]])
        >>> df = DepressionFinderAndRouter(rg, routing="D8")
        >>> df.map_depressions()
        >>> rg.at_node["flow__receiver_node"].reshape(rg.shape)
        array([[ 0,  1,  2,  3,  4,  5,  6],
               [ 7,  7, 16, 16, 10, 18, 13],
               [14, 14,  8, 16, 17, 18, 20],
               [21, 21, 16, 16, 24, 33, 27],
               [28, 28, 24, 24, 24, 32, 34],
               [35, 35, 38, 32, 38, 32, 41],
               [42, 43, 44, 45, 46, 47, 48]])
        """

        (links, nbrs, diag_nbrs) = self._links_and_nbrs_at_node(outlet_node)

        # Sweep through them, identifying the neighbor with the greatest slope.
        # We are probably duplicating some gradient calculations, but this only
        # happens occasionally, when we have a candidate outlet node.

        # We're seeking the valid neighbor (*receiver*) in the direction of
        # steepest descent. Initially set *receiver* to the node itself, and
        # downhill-positive gradient to zero. If we don't find any neighbor
        # with a steeper path (or an open boundary), then we have failed.
        max_downhill_grad = 0.0
        receiver = outlet_node
        node_elev = self._elev[outlet_node]

        # Iterate over all "regular" neighbors
        for i in range(len(links)):
            lnk = links[i]
            nbr = nbrs[i]

            # To pass this first hurdle, the neighbor must:
            #   * not be part of the current lake
            #   * have a surface (if flooded, WATER surface)
            #     lower than our outlet node;
            #   * not be a closed boundary
            if (
                self._flood_status[nbr] != _CURRENT_LAKE
                and (
                    (self._elev[nbr] + self._depression_depth[nbr])
                    < self._elev[receiver]
                )
                and self._grid.status_at_node[nbr] != self._grid.BC_NODE_IS_CLOSED
            ):
                # Next test: is it the steepest downhill grad so far?
                # If so, we've found a candidate.
                grad = (node_elev - self._elev[nbr]) / self._grid.length_of_link[lnk]
                if grad > max_downhill_grad:
                    # Update the receiver and max grad: this is now the one
                    # to beat.
                    max_downhill_grad = grad
                    receiver = nbr

        # If we're on a D8 raster, iterate over all diagonal neighbors
        if self._D8:
            for nbr in diag_nbrs:
                # Again, to pass this first hurdle, the neighbor must:
                #   * not be part of the current lake
                #   * have a surface (if flooded, WATER surface)
                #     lower than our outlet node;
                #   * not be a closed boundary
                if (
                    self._flood_status[nbr] != _CURRENT_LAKE
                    and (
                        (self._elev[nbr] + self._depression_depth[nbr])
                        < self._elev[receiver]
                    )
                    and self._grid.status_at_node[nbr] != self._grid.BC_NODE_IS_CLOSED
                ):
                    # Next test: is it the steepest downhill grad so far?
                    # If so, we've found a candidate.
                    grad = (node_elev - self._elev[nbr]) / self._diag_link_length
                    if grad > max_downhill_grad:
                        # Update the receiver and max grad: this is now the one
                        # to beat.
                        max_downhill_grad = grad
                        receiver = nbr

        # We only call this method after is_valid_outlet has evaluated True,
        # so in theory it should NEVER be the case that we fail to find a
        # receiver. However, let's make sure.
        assert receiver != outlet_node, "failed to find receiver with ID: %r" % receiver

        # Finally, let's assign it

        self._grid.at_node["flow__receiver_node"][outlet_node] = receiver

    def node_can_drain(self, the_node):
        """Check if a node has drainage away from the current lake/depression.

        Parameters
        ----------
        the_node : int
            The node to test.
        nodes_this_depression : array_like of int
            Nodes that form a pit.

        Returns
        -------
        boolean
            ``True`` if the node can drain. Otherwise, ``False``.
        """
        nbrs = self._node_nbrs[the_node]
        not_bad = nbrs != self._grid.BAD_INDEX
        not_too_high = self._elev[nbrs] < self._elev[the_node]
        not_current_lake = np.not_equal(self._flood_status[nbrs], _CURRENT_LAKE)
        not_flooded = np.not_equal(self._flood_status[nbrs], _FLOODED)

        # The following logic block handles the case when a neighbor is
        # flooded but its outlet is LOWER than the_node, so the_node could
        # be an outlet that flows into a lower lake.
        #
        # We proceed only if there is at least one flooded node
        if np.any(np.logical_not(not_flooded)):
            # Examine each neighbor
            for i in range(len(nbrs)):
                # If the neighbor is flooded...
                if not not_flooded[i]:
                    # Check to see whether its own outlet is lower than
                    # the_node. If so, then it does not "count" as being
                    # flooded, because its water level is lower than our
                    # current potential lake outlet.
                    dep_out = self._depression_outlet_map[nbrs[i]]
                    if self._elev[the_node] > self._elev[dep_out]:
                        not_flooded[i] = True

        # Now combine all the issues: any neighbor(s) that is not "bad",
        # too high, part of the current lake, or flooded at a level equal to
        # or higher than the_node, is a potential outlet. So, if there are any
        # neighbor nodes that pass all these tests, then the_node can drain.
        all_probs = np.logical_and(
            np.logical_and(not_bad, not_too_high),
            np.logical_and(not_current_lake, not_flooded),
        )
        return np.any(all_probs)

    def is_valid_outlet(self, the_node):
        """Check if a node is a valid outlet for the depression.

        Parameters
        ----------
        the_node : int
            The node to test.
        nodes_this_depression : array_like of int
            Nodes that form a pit.

        Returns
        -------
        boolean
            ``True`` if the node is a valid outlet. Otherwise, ``False``.
        """
        if self._grid.status_at_node[the_node] == self._grid.BC_NODE_IS_FIXED_VALUE:
            return True

        if self.node_can_drain(the_node):
            return True

        return False

    def _record_depression_depth_and_outlet(
        self, nodes_this_depression, outlet_id, pit_node
    ):
        """Record information about a depression.

        Record information about this depression/lake in the flood_status,
        depression_depth, and depression_outlet arrays.

        Parameters
        ----------
        nodes_this_depression : iterable of int
            Nodes that form a pit.
        outlet_id : int
            Node that is the outlet of the pit.
        pit_node : int
            Node that is the deepest pit, uniquely associated with this
            depression.
        """
        n = nodes_this_depression

        # three cases possible - new lake is fresh; new lake is smaller than
        # an existing lake (subsumed, and unimportant), new lake is equal to
        # or bigger than old lake (or multiple old lakes). It SHOULDN'T be
        # possible to have two lakes overlapping... We can test this with an
        # assertion that out total # of *tracked* lakes matches the accumulated
        # total of unique vals in lake_map.
        fresh_nodes = np.equal(self._lake_map[n], self._grid.BAD_INDEX)
        if np.all(fresh_nodes):  # a new lake
            self._flood_status[n] = _FLOODED
            self._depression_depth[n] = self._elev[outlet_id] - self._elev[n]
            self._depression_outlet_map[n] = outlet_id
            self._lake_map[n] = pit_node
            self._pits_flooded += 1
            pit_node_where = np.searchsorted(self._pit_node_ids, pit_node)
            self._unique_pits[pit_node_where] = True
        elif np.any(fresh_nodes):  # lake is bigger than one or more existing
            self._flood_status[n] = _FLOODED
            depth_this_lake = self._elev[outlet_id] - self._elev[n]
            self._depression_depth[n] = depth_this_lake
            self._depression_outlet_map[n] = outlet_id
            # ^these two will just get stamped over as needed
            subsumed_lakes = np.unique(self._lake_map[n])  # IDed by pit_node
            # the final entry is self._grid.BAD_INDEX
            subs_lakes_where = np.searchsorted(self._pit_node_ids, subsumed_lakes[1:])
            pit_node_where = np.searchsorted(self._pit_node_ids, pit_node)
            self._unique_pits[subs_lakes_where] = False
            self._unique_pits[pit_node_where] = True
            self._pits_flooded -= subsumed_lakes.size - 2
            # -1 for the self._grid.BAD_INDEX that must be present; another -1
            # because a single lake is just replaced by a new lake.
            self._lake_map[n] = pit_node
        else:  # lake is subsumed within an existing lake
            print(" eaten lake")
            assert np.all(np.equal(self._flood_status[n], _CURRENT_LAKE))
            self._flood_status[n] = _FLOODED

    def find_depression_from_pit(self, pit_node, reroute_flow=True):
        """Find the extent of the nodes that form a pit.

        Identify extent of depression/lake whose lowest point is the node
        pit_node (which is a itself a pit, a.k.a., closed depression).

        Parameters
        ----------
        pit_node : int
            The node that is the lowest point of a pit.
        """
        # Flag the pit as being _CURRENT_LAKE (it's the first node in the
        # current lake)
        self._flood_status[pit_node] = _CURRENT_LAKE

        # This flag keeps track of when we're done with this depression
        found_outlet = False

        # Safety check
        count = 0
        max_count = self._grid.number_of_nodes + 1

        # Place pit_node at top of depression list
        nodes_this_depression = self.grid.zeros("node", dtype=int)
        nodes_this_depression[0] = pit_node
        pit_count = 1

        while not found_outlet:
            lowest_node_on_perimeter, pit_count = find_lowest_node_on_lake_perimeter_c(
                self._node_nbrs,
                self.flood_status,
                self._elev,
                nodes_this_depression,
                pit_count,
                self._BIG_ELEV,
            )
            # note this can return the supplied node, if - somehow - the
            # surrounding nodes are all self._grid.BAD_INDEX
            # I BELIEVE THE IS_VALID_OUTLET FN SHOULD ASSIGN FLOW DIR
            found_outlet = self.is_valid_outlet(lowest_node_on_perimeter)

            # If we haven't found an outlet, add lowest_node to the lake list
            # and flag it as being part of the current lake/depression
            if not found_outlet:
                nodes_this_depression[pit_count] = lowest_node_on_perimeter
                self._flood_status[lowest_node_on_perimeter] = _CURRENT_LAKE
                pit_count += 1

            # If we HAVE found an outlet, and we are re-routing flow, then
            # assign the proper flow direction to the outlet node. If it is an
            # open boundary, then it drains to itself. Otherwise, call
            # assign_outlet_receiver to find the correct receiver (so that it
            # doesn't simply drain back into the lake)
            elif ("flow__receiver_node" in self._grid.at_node) and reroute_flow:
                if (
                    self._grid.status_at_node[lowest_node_on_perimeter]
                    != self._grid.BC_NODE_IS_CORE
                ):
                    self._grid.at_node["flow__receiver_node"][
                        lowest_node_on_perimeter
                    ] = lowest_node_on_perimeter
                else:
                    self.assign_outlet_receiver(lowest_node_on_perimeter)

            # Safety check, in case a bug (ha!) puts us in an infinite loop
            assert count < max_count, "too many iterations in lake filler!"
            count += 1

        self._depression_outlets.append(lowest_node_on_perimeter)
        # Now that we've mapped this depression, record it in the arrays
        # depression_depth, depression_outlet, and flood_status
        self._record_depression_depth_and_outlet(
            nodes_this_depression[:pit_count], lowest_node_on_perimeter, pit_node
        )

        # TODO: ideally we need a way to keep track of the number, area extent,
        # and average depth of depressions. Tricky thing is that one might be
        # devoured by another, so would need to be removed from the list.

    def _identify_depressions_and_outlets(self, reroute_flow=True):
        """Find depression and lakes on a topographic surface.

        Find and map the depressions/lakes in a topographic surface,
        given a previously identified list of pits (if any) in the
        surface.
        """
        self._pits_flooded = 0
        self._unique_pits = np.zeros_like(self._pit_node_ids, dtype=bool)
        # debug_count = 0
        for pit_node in self._pit_node_ids:
            if self._flood_status[pit_node] != _PIT:
                self._depression_outlets.append(self._grid.BAD_INDEX)
            else:
                self.find_depression_from_pit(pit_node, reroute_flow)
                self._pits_flooded += 1

        assert len(self._depression_outlets) == self._unique_pits.size

        self._unique_lake_outlets = np.array(self._depression_outlets)[
            self._unique_pits
        ]

    def update(self):
        """Alias for map_depressions."""
        self.map_depressions()

    def map_depressions(self):
        """Map depressions/lakes in a topographic surface.

        Examples
        --------
        Test #1: 5x5 raster grid with a diagonal lake.

        >>> import numpy as np
        >>> from landlab import RasterModelGrid
        >>> from landlab.components import DepressionFinderAndRouter

        >>> rg = RasterModelGrid((5, 5))
        >>> rg.at_node["topographic__elevation"] = [
        ...     [100.0, 100.0, 95.0, 100.0, 100.0],
        ...     [100.0, 101.0, 92.0, 1.0, 100.0],
        ...     [100.0, 101.0, 2.0, 101.0, 100.0],
        ...     [100.0, 3.0, 101.0, 101.0, 100.0],
        ...     [90.0, 95.0, 100.0, 100.0, 100.0],
        ... ]

        >>> df = DepressionFinderAndRouter(rg, reroute_flow=False)
        >>> df.map_depressions()
        >>> df.display_depression_map()
        . . . . .
        . . . ~ .
        . . ~ . .
        . ~ . . .
        o . . . .
        """
        if self._bc_set_code != self._grid.bc_set_code:
            self.updated_boundary_conditions()
            self._bc_set_code = self._grid.bc_set_code

        # verify that there is an outlet to the grid and
        if not np.any(
            self._grid.status_at_node[self._grid.boundary_nodes]
            != self._grid.BC_NODE_IS_CLOSED
        ):
            raise ValueError(
                "DepressionFinderAndRouter requires that there is at least one "
                "open boundary node."
            )

        self._lake_map.fill(self._grid.BAD_INDEX)
        self._depression_outlet_map.fill(self._grid.BAD_INDEX)
        self._depression_depth.fill(0.0)
        self._depression_outlets = []  # reset these
        # Locate nodes with pits
        if isinstance(self._user_supplied_pits, str):
            try:
                pits = self._grid.at_node[self._user_supplied_pits]
                supplied_pits = np.where(pits)[0]
                self._pit_node_ids = as_id_array(
                    np.setdiff1d(supplied_pits, self._grid.boundary_nodes)
                )
                self._number_of_pits = self._pit_node_ids.size
                self._is_pit.fill(False)
                self._is_pit[self._pit_node_ids] = True
            except FieldError:
                self._find_pits()
        elif self._user_supplied_pits is None:
            self._find_pits()
        else:  # hopefully an array or other sensible iterable
            if len(self._user_supplied_pits) == self._grid.number_of_nodes:
                supplied_pits = np.where(self._user_supplied_pits)[0]
            else:  # it's an array of node ids
                supplied_pits = self._user_supplied_pits
            # remove any boundary nodes from the supplied pit list
            self._pit_node_ids = as_id_array(
                np.setdiff1d(supplied_pits, self._grid.boundary_nodes)
            )

            self._number_of_pits = self._pit_node_ids.size
            self._is_pit.fill(False)
            self._is_pit[self._pit_node_ids] = True
        # Set up "lake code" array
        self._flood_status.fill(_UNFLOODED)
        self._flood_status[self._pit_node_ids] = _PIT

        self._identify_depressions_and_outlets(self._reroute_flow)

        if self._reroute_flow and ("flow__receiver_node" in self._grid.at_node):
            self._receivers = self._grid.at_node["flow__receiver_node"]
            self._sinks = self._grid.at_node["flow__sink_flag"]
            self._grads = self._grid.at_node["topographic__steepest_slope"]
            self._links = self._grid.at_node["flow__link_to_receiver_node"]
            self._route_flow()
            self._reaccumulate_flow()

    def _find_unresolved_neighbors(self, nbrs, receivers):
        """Make and return list of neighbors of node with unresolved flow dir.

        Examples
        --------
        >>> import numpy as np
        >>> from landlab.components import DepressionFinderAndRouter
        >>> from landlab import RasterModelGrid
        >>> rg = RasterModelGrid((7, 8))
        >>> z = rg.add_zeros("topographic__elevation", at="node")
        >>> df = DepressionFinderAndRouter(rg)
        >>> rcvr = np.arange(56)
        >>> rcvr[13] = -1
        >>> rcvr[21] = -1
        >>> rcvr[29] = -1
        >>> rcvr[30] = -1
        >>> nbrs = np.array([23, 30, 21, 14], dtype=int)
        >>> df._find_unresolved_neighbors(nbrs, rcvr)
        array([30, 21])
        """
        # unresolved = np.where(receivers[nbrs] == -1)[0]
        # ur_nbrs = nbrs[unresolved]
        # ur_links = self._grid.links_at_node[unresolved]
        # return (ur_nbrs, ur_links)
        return nbrs[np.where(receivers[nbrs] == -1)[0]]

    def _find_unresolved_neighbors_new(self, nbrs, nbr_links, receivers):
        """Make and return list of neighbors of node with unresolved flow dir.

        Examples
        --------
        >>> import numpy as np
        >>> from landlab.components import DepressionFinderAndRouter
        >>> from landlab import RasterModelGrid
        >>> rg = RasterModelGrid((7, 8))
        >>> z = rg.add_zeros("topographic__elevation", at="node")
        >>> df = DepressionFinderAndRouter(rg)
        >>> rcvr = np.arange(56)
        >>> rcvr[13] = -1
        >>> rcvr[21] = -1
        >>> rcvr[29] = -1
        >>> rcvr[30] = -1
        >>> nbrs = rg.adjacent_nodes_at_node[22]
        >>> nbr_links = rg.links_at_node[22]
        >>> df._find_unresolved_neighbors_new(nbrs, nbr_links, rcvr)
        (array([30, 21]), array([43, 35]))
        >>> nbrs = rg.diagonal_adjacent_nodes_at_node[22]
        >>> nbr_links = rg.d8s_at_node[22, 4:]
        >>> df._find_unresolved_neighbors_new(nbrs, nbr_links, rcvr)
        (array([29, 13]), array([136, 121]))
        """
        unresolved = np.where(receivers[nbrs] == -1)[0]
        ur_nbrs = nbrs[unresolved]
        ur_links = nbr_links[unresolved]
        return (ur_nbrs, ur_links)

    def _route_flow_for_one_lake(self, outlet, lake_nodes):
        """Route flow across a single lake. Alternative to part of _route_flow.

        Examples
        --------
        >>> from landlab import RasterModelGrid
        >>> import numpy as np
        >>> from landlab.components import DepressionFinderAndRouter
        >>> rg = RasterModelGrid((7, 8))
        >>> z = rg.add_zeros("topographic__elevation", at="node")
        >>> rcvr = rg.add_zeros("flow__receiver_node", at="node", dtype=int)
        >>> rcvr[:] = np.arange(rg.number_of_nodes)
        >>> lake_nodes = np.flatnonzero(
        ...     [
        ...         [0, 0, 0, 0, 0, 0, 0, 0],
        ...         [0, 0, 1, 0, 1, 1, 0, 0],
        ...         [0, 0, 0, 1, 1, 1, 0, 0],
        ...         [0, 1, 1, 1, 1, 1, 1, 0],
        ...         [0, 1, 1, 1, 1, 1, 1, 0],
        ...         [0, 0, 0, 0, 1, 1, 1, 0],
        ...         [0, 0, 0, 0, 0, 0, 0, 0],
        ...     ]
        ... )

        >>> rcvr[9] = 1
        >>> rcvr[11] = 3
        >>> rcvr[14] = 6
        >>> rcvr[17] = 16
        >>> rcvr[18] = 17
        >>> rcvr[22] = 14  # this is the outlet
        >>> rcvr[41] = 40
        >>> rcvr[42] = 50
        >>> rcvr[43] = 51
        >>> df = DepressionFinderAndRouter(rg)
        >>> df.receivers = rcvr
        >>> df._route_flow_for_one_lake(22, lake_nodes)
        >>> df.receivers
        array([ 0,  1,  2,  3,  4,  5,  6,  7,  8,  1, 19,  3, 13, 22,  6, 15, 16,
               16, 17, 20, 21, 22, 14, 23, 24, 26, 27, 28, 29, 22, 22, 31, 32, 34,
               35, 36, 29, 29, 30, 39, 40, 40, 50, 51, 36, 37, 38, 47, 48, 49, 50,
               51, 52, 53, 54, 55])
        """

        # Flag receiver nodes inside the lake as "unresolved"
        UNRESOLVED = -1
        self._receivers[lake_nodes] = UNRESOLVED

        # We work with two lists: the nodes currently being processed, and
        # the nodes that will be processed on the next iteration. We start with
        # the outlet node as the one being processed, and an empty list of
        # nodes to process next.
        nodes_being_processed = [outlet]
        nodes_to_proc_next = []

        # We must now iterate until we've taken care of all the nodes in the
        # lake. In each iteration, we:
        #  1 - find the unresolved neighbors of nodes being processed
        #  2 - point them toward the nodes being processed
        #  3 - place them on the nodes_to_proc_next list
        # We stop when there are no more nodes to process.
        #    Note that the nested looping will be slow, but could be sped up
        # by translating to cython.
        counter = 0  # counts # of times thru loop as fail-safe
        done = False
        while not done:
            # Get unresolved "regular" neighbors of the current nodes
            for cn in nodes_being_processed:
                # Get active and unresolved neighbors of cn
                (nbrs, lnks) = self._find_unresolved_neighbors_new(
                    self._grid.adjacent_nodes_at_node[cn],
                    self._grid.links_at_node[cn],
                    self._receivers,
                )
                # They will now flow to cn
                if nbrs.size > 0:
                    self._receivers[nbrs] = cn
                    if "flow__link_to_receiver_node" in self._grid.at_node:
                        self._links[nbrs] = lnks
                        slopes = (
                            self._elev[nbrs] - self._elev[cn]
                        ) / self._grid.length_of_link[lnks]
                        self._grads[nbrs] = np.maximum(slopes, 0.0)

                # Place them on the list of nodes to process next
                for n in nbrs:
                    nodes_to_proc_next.append(n)

            # If we're working with a raster that has diagonals, do the same
            # for the diagonal neighbors
            if self._D8:
                # Get unresolved "regular" neighbors of the current nodes
                for cn in nodes_being_processed:
                    (nbrs, diags) = self._find_unresolved_neighbors_new(
                        self._grid.diagonal_adjacent_nodes_at_node[cn],
                        self._grid.d8s_at_node[cn, 4:],
                        self._receivers,
                    )
                    # They will now flow to cn
                    if nbrs.size > 0:
                        self._receivers[nbrs] = cn
                        if "flow__link_to_receiver_node" in self._grid.at_node:
                            self._links[nbrs] = diags
                            slopes = (
                                self._elev[nbrs] - self._elev[cn]
                            ) / self._diag_link_length
                            self._grads[nbrs] = np.maximum(slopes, 0.0)

                    # Place them on the list of nodes to process next
                    for n in nbrs:
                        nodes_to_proc_next.append(n)

            # Move to the next set of nodes
            nodes_being_processed = nodes_to_proc_next
            nodes_to_proc_next = []
            if not nodes_being_processed:
                done = True

            # Just in case
            counter += 1
            assert counter < self._grid.number_of_nodes, "inf loop in lake"

    def _route_flow(self):
        """Route flow across lake flats.

        Route flow across lake flats, which have already been
        identified.
        """

        # Process each lake.
        for outlet_node, lake_code in zip(self.lake_outlets, self.lake_codes):
            # Get the nodes in the lake
            nodes_in_lake = np.where(self._lake_map == lake_code)[0]
            if len(nodes_in_lake) > 0:
                # find the correct outlet for the lake, if necessary
                if self._lake_map[self._receivers[outlet_node]] == lake_code:
                    nbrs = self._grid.active_adjacent_nodes_at_node[outlet_node]
                    not_lake = nbrs[np.where(self._lake_map[nbrs] != lake_code)[0]]
                    min_index = np.argmin(self._elev[not_lake])
                    new_receiver = not_lake[min_index]

                    # set receiver for new outlet.
                    self._receivers[outlet_node] = new_receiver

                # reset_link for new outlet
                outlet_receiver = self._receivers[outlet_node]
                if self._D8:
                    adjacent_links_and_diags = np.hstack(
                        (
                            self._grid.adjacent_nodes_at_node[outlet_node, :],
                            self._grid.diagonal_adjacent_nodes_at_node[outlet_node, :],
                        )
                    )
                    find_recs = outlet_receiver == adjacent_links_and_diags
                    new_link = self._grid.d8s_at_node[outlet_node, find_recs]
                else:
                    find_recs = (
                        outlet_receiver
                        == self._grid.adjacent_nodes_at_node[outlet_node, :]
                    )
                    new_link = self._grid.links_at_node[outlet_node, find_recs]

                if new_link.size == 0:
                    new_link = self._grid.BAD_INDEX
                if np.min(new_link) == np.max(new_link) and np.min(new_link) == -1:
                    self._links[outlet_node] = -1
                else:
                    self._links[outlet_node] = new_link

                # make a check
                assert (
                    self._lake_map[self._receivers[outlet_node]] != lake_code
                ), "outlet of lake drains to itself!"

                # Route flow
                self._route_flow_for_one_lake(outlet_node, nodes_in_lake)

        self._sinks[self._pit_node_ids] = False

    def _reaccumulate_flow(self):
        """Update drainage area, discharge, upstream order, and flow link.

        Invoke the accumulator a second time to update drainage area,
        discharge, and upstream order.
        """
        # Calculate drainage area, discharge, and downstr->upstr order
        Q_in = self._grid.at_node["water__unit_flux_in"]
        areas = self._grid.cell_area_at_node.copy()
        areas[self._grid.closed_boundary_nodes] = 0.0

        self._a, q, s = flow_accum_bw.flow_accumulation(
            self._receivers, node_cell_area=areas, runoff_rate=Q_in
        )

        # finish the property updating:
        self._grid.at_node["drainage_area"][:] = self._a
        self._grid.at_node["surface_water__discharge"][:] = q
        self._grid.at_node["flow__upstream_node_order"][:] = s

    def _handle_outlet_node(self, outlet_node, nodes_in_lake):
        """Ensure the outlet node drains to the grid edge.

        Makes sure the outlet node is drains to the grid edge, not back
        into the depression.
        This exists because if the slope into the lake is steeper than the
        slope out from the (rim lowest) outlet node, the lake won't drain.

        Parameters
        ----------
        outlet_node : int
            The outlet node.
        nodes_in_lake : array_like of int
            The nodes that are contained within the lake.
        """
        if self._grid.status_at_node[outlet_node] == 0:  # it's not a BC
            if self._D8:
                outlet_neighbors = np.hstack(
                    (
                        self._grid.active_adjacent_nodes_at_node[outlet_node],
                        self._grid.diagonal_adjacent_nodes_at_node[outlet_node],
                    )
                )
            else:
                outlet_neighbors = self._grid.active_adjacent_nodes_at_node[
                    outlet_node
                ].copy()
            inlake = np.in1d(outlet_neighbors.flat, nodes_in_lake)
            assert inlake.size > 0
            outlet_neighbors[inlake] = -1
            unique_outs, unique_indxs = np.unique(outlet_neighbors, return_index=True)
            out_draining = unique_outs[1:]
            if isinstance(self._grid, RasterModelGrid):
                link_l = self._link_lengths
            else:  # Voronoi
                link_l = self._link_lengths[self._grid.links_at_node[outlet_node, :]]
            eff_slopes = (self._elev[outlet_node] - self._elev[out_draining]) / link_l[
                unique_indxs[1:]
            ]
            lowest = np.argmax(eff_slopes)
            lowest_node = out_draining[lowest]
            # route the flow
            self._receivers[outlet_node] = lowest_node
        else:
            self._receivers[outlet_node] = outlet_node

    def display_depression_map(self):
        """Print a simple character-based map of depressions/lakes."""
        # Find the outlet nodes (just for display purposes)
        is_outlet = np.zeros(self._grid.number_of_nodes, dtype=bool)
        for i in self._grid.core_nodes:
            if self._flood_status[i] == _FLOODED:
                is_outlet[self._depression_outlet_map[i]] = True

        n = 0
        for _ in range(self._grid.number_of_node_rows):
            for _ in range(self._grid.number_of_node_columns):
                if is_outlet[n]:
                    print("o", end=" ")
                elif self._flood_status[n] == _UNFLOODED:
                    print(".", end=" ")
                else:
                    print("~", end=" ")
                n += 1
            print()

    @property
    def lake_outlets(self):
        """Returns the *unique* outlets for each lake, in same order as the
        return from lake_codes."""
        return np.array(self._depression_outlets)[self._unique_pits]

    @property
    def lake_codes(self):
        """Returns the *unique* code assigned to each unique lake.

        These are the values used to map the lakes in the property
        "lake_map".
        """
        return self._pit_node_ids[self._unique_pits]

    @property
    def number_of_lakes(self):
        """Return the number of individual lakes."""
        return self._unique_pits.sum()

    @property
    def lake_map(self):
        """Return an array of ints, where each node within a lake is labelled
        with a unique (non-consecutive) code corresponding to each unique lake.

        The codes used can be obtained with *lake_codes*. Nodes not in a
        lake are labelled with self._grid.BAD_INDEX
        """
        return self._lake_map

    @property
    def lake_at_node(self):
        """Return a boolean array, True if the node is flooded, False
        otherwise."""
        return self._lake_map != self._grid.BAD_INDEX

    @property
    def lake_areas(self):
        """A nlakes-long array of the area of each lake.

        The order is the same as that returned by *lake_codes*.
        """
        lake_areas = np.empty(self.number_of_lakes)
        for lake_counter, lake_code in enumerate(self.lake_codes):
            each_cell_in_lake = self._grid.cell_area_at_node[
                self._lake_map == lake_code
            ]
            lake_areas[lake_counter] = each_cell_in_lake.sum()
        return lake_areas

    @property
    def lake_volumes(self):
        """A nlakes-long array of the volume of each lake.

        The order is the same as that returned by *lake_codes*.
        """
        lake_vols = np.empty(self.number_of_lakes)
        col_vols = self._grid.cell_area_at_node * self._depression_depth
        for lake_counter, lake_code in enumerate(self.lake_codes):
            each_cell_in_lake = col_vols[self._lake_map == lake_code]
            lake_vols[lake_counter] = each_cell_in_lake.sum()
        return lake_vols

"""Find depressions on a topographic surface.

.. codeauthor:: gtucker, DEJH (Flow routing)
"""

# Routing by DEJH, Oct 15.
import warnings
from io import StringIO

import numpy as np

from landlab.grid.nodestatus import NodeStatus

from ...core.model_component import Component
from ...core.utils import as_id_array
from ...field import FieldError
from ...grid import RasterModelGrid
from ..flow_accum import flow_accum_bw
from .cfuncs import find_lowest_node_on_lake_perimeter_c
from .cfuncs import find_pits
from .errors import NoOutletError
from .floodstatus import FloodStatus


class DepressionFinderAndRouter(Component):
    """Find depressions on a topographic surface.

    This component identifies depressions in a topographic surface, finds an
    outlet for each depression [1]_.  If directed to do so (the default), and the
    component is able to find existing routing fields created by the
    `route_flow_dn` component, it will then modify the drainage directions and
    accumulations already stored in the grid to route flow across these
    depressions.

    Note that, in general, properties of this class named *"depression"* identify
    each individual pit in the topography, including those that will merge
    once the fill is performed. Those named *"lake"* return the unique lakes
    created by the fill, and are probably the properties most users will
    want.

    Note also that the structure of drainage within the lakes is not
    guaranteed, and in particular, may not be symmetrical even if your
    boundary conditions are. However, the outputs from the lake will all
    still be correct.

    Note the routing part of this component may not yet be fully compatible with
    irregular grids.

    The prinary method of this class is `map_depressions`.

    Examples
    --------
    Route flow across a depression in a sloped surface.

    >>> from landlab import RasterModelGrid
    >>> from landlab.components import FlowAccumulator, DepressionFinderAndRouter
    >>> from landlab.grid.raster_funcs import neighbor_to_arrow

    >>> grid = RasterModelGrid((7, 7), xy_spacing=1.0)
    >>> z = (grid.x_of_node * 100.0 + grid.y_of_node).reshape(grid.shape)
    >>> z[2:-2, 2:-2] *= 0.1
    >>> z
    array([[  0. , 100. , 200. , 300. , 400. , 500. , 600. ],
           [  1. , 101. , 201. , 301. , 401. , 501. , 601. ],
           [  2. , 102. ,  20.2,  30.2,  40.2, 502. , 602. ],
           [  3. , 103. ,  20.3,  30.3,  40.3, 503. , 603. ],
           [  4. , 104. ,  20.4,  30.4,  40.4, 504. , 604. ],
           [  5. , 105. , 205. , 305. , 405. , 505. , 605. ],
           [  6. , 106. , 206. , 306. , 406. , 506. , 606. ]])
    >>> grid.at_node["topographic__elevation"] = z

    >>> fr = FlowAccumulator(grid, flow_director="D8")
    >>> fr.run_one_step()
    >>> grid.at_node["flow__receiver_node"].reshape(grid.shape)
    array([[ 0,  1,  2,  3,  4,  5,  6],
           [ 7,  7, 16, 17, 18, 18, 13],
           [14, 14, 16, 16, 17, 18, 20],
           [21, 21, 16, 23, 24, 25, 27],
           [28, 28, 23, 30, 31, 32, 34],
           [35, 35, 30, 31, 32, 32, 41],
           [42, 43, 44, 45, 46, 47, 48]])

    The flow gets stuck in a depression at (2, 2).

    >>> neighbor_to_arrow(grid.at_node["flow__receiver_node"].reshape(grid.shape))
    array([['○', '○', '○', '○', '○', '○', '○'],
           ['○', '←', '↓', '↓', '↓', '↙', '○'],
           ['○', '←', '○', '←', '←', '←', '○'],
           ['○', '←', '↑', '←', '←', '←', '○'],
           ['○', '←', '↑', '←', '←', '←', '○'],
           ['○', '←', '↑', '↑', '↑', '↖', '○'],
           ['○', '○', '○', '○', '○', '○', '○']], dtype='<U1')

    >>> grid.at_node["drainage_area"].reshape(grid.shape)
    array([[ 0.,  0.,  0.,  0.,  0.,  0.,  0.],
           [ 1.,  1.,  1.,  1.,  1.,  1.,  0.],
           [ 1.,  1., 20.,  6.,  4.,  1.,  0.],
           [ 1.,  1., 12.,  3.,  2.,  1.,  0.],
           [ 1.,  1.,  8.,  6.,  4.,  1.,  0.],
           [ 1.,  1.,  1.,  1.,  1.,  1.,  0.],
           [ 0.,  0.,  0.,  0.,  0.,  0.,  0.]])

    >>> df = DepressionFinderAndRouter(grid, routing="D8", reroute_flow=True)
    >>> df.map_depressions()

    The flow now is no longer stuck and the lake now drains to node 8.

    >>> df.lake_outlets
    array([8])
    >>> neighbor_to_arrow(grid.at_node["flow__receiver_node"].reshape(grid.shape))
    array([['○', '○', '○', '○', '○', '○', '○'],
           ['○', '←', '↓', '↓', '↓', '↙', '○'],
           ['○', '←', '↖', '←', '←', '←', '○'],
           ['○', '←', '↑', '↖', '←', '←', '○'],
           ['○', '←', '↑', '↑', '↖', '←', '○'],
           ['○', '←', '↑', '↑', '↑', '↖', '○'],
           ['○', '○', '○', '○', '○', '○', '○']], dtype='<U1')


    The area of each lake in lake.

    >>> df.lake_areas
    array([9.])
    >>> df.lake_at_node.reshape(grid.shape).view(dtype=np.uint8)
    array([[0, 0, 0, 0, 0, 0, 0],
           [0, 0, 0, 0, 0, 0, 0],
           [0, 0, 1, 1, 1, 0, 0],
           [0, 0, 1, 1, 1, 0, 0],
           [0, 0, 1, 1, 1, 0, 0],
           [0, 0, 0, 0, 0, 0, 0],
           [0, 0, 0, 0, 0, 0, 0]], dtype=uint8)
    >>> grid.at_node["drainage_area"].reshape(grid.shape)
    array([[ 0.,  0.,  0.,  0.,  0.,  0.,  0.],
           [21., 21.,  1.,  1.,  1.,  1.,  0.],
           [ 1.,  1., 20.,  6.,  4.,  1.,  0.],
           [ 1.,  1.,  3.,  9.,  2.,  1.,  0.],
           [ 1.,  1.,  2.,  2.,  4.,  1.,  0.],
           [ 1.,  1.,  1.,  1.,  1.,  1.,  0.],
           [ 0.,  0.,  0.,  0.,  0.,  0.,  0.]])

    A unique code for each lake present on the grid.

    >>> df.lake_codes
    array([16])
    >>> df.lake_map.reshape(grid.shape)
    array([[-1, -1, -1, -1, -1, -1, -1],
           [-1, -1, -1, -1, -1, -1, -1],
           [-1, -1, 16, 16, 16, -1, -1],
           [-1, -1, 16, 16, 16, -1, -1],
           [-1, -1, 16, 16, 16, -1, -1],
           [-1, -1, -1, -1, -1, -1, -1],
           [-1, -1, -1, -1, -1, -1, -1]])

    Because `reroute_flow` was set to ``True``, the flow connectivity fields
    created by the :py:class:`~landlab.components.flow_accum.FlowAccumulator`
    will have now been modified to route flow over the depressions in the
    surface. The topogrphy itself is not modified.

    References
    ----------
    **Required Software Citation(s) Specific to this Component**

    None Listed

    **Additional References**

    .. [1] Tucker, G. E., Lancaster, S. T., Gasparini, N. M., and Bras, R. L.: The
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
        "topographic__elevation": {
            "dtype": float,
            "intent": "in",
            "optional": False,
            "units": "m",
            "mapping": "node",
            "doc": "Land surface topographic elevation",
        },
    }

    def __init__(self, grid, routing=None, pits="flow__sink_flag", reroute_flow=True):
        """Create a DepressionFinderAndRouter.

        Constructor assigns a copy of the grid, sets the current time, and
        calls the initialize method.

        Parameters
        ----------
        grid : RasterModelGrid
            A landlab RasterModelGrid.
        routing : {'D4', 'D8'}, optional
            If grid is a raster type, controls whether lake connectivity can
            occur on diagonals ('D8', default), or only orthogonally ('D4').
            Has no effect if grid is not a raster.
        pits : array_like or str, optional
            If a field name, the boolean field containing True where pits.
            If array-like, either a boolean array of nodes of the pits, or an
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

        if pits is None or isinstance(pits, str):
            self._pits = pits
        else:
            self._pits = np.asarray(pits)
        self._reroute_flow = reroute_flow

        self._routing = None
        if isinstance(grid, RasterModelGrid):
            self._routing = self._validate_routing(routing) if routing else "D8"
        elif routing in {"D4", "D8"}:
            warnings.warn(
                f"ignoring supplied routing method {routing!r}, not a raster grid",
                stacklevel=2,
            )

        if "flow__receiver_node" in grid.at_node and grid.at_node[
            "flow__receiver_node"
        ].size != grid.size("node"):
            raise NotImplementedError(
                "A route-to-multiple flow director has been"
                " run on this grid. The depression finder is"
                " not compatible with the grid anymore. Use"
                " DepressionFinderAndRouter with reroute_flow=True"
                " only with route-to-one methods. If using this"
                " component with such a flow directing method is desired,"
                " please open a GitHub Issue"
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

        self.updated_boundary_conditions()

        self._lake_outlets = []  # a list of each unique lake outlet
        # ^note this is nlakes-long

        self._flood_status = np.zeros(self.grid.number_of_nodes, dtype=int)
        self._pit_nodes = None

        self._lake_map = np.empty(self._grid.number_of_nodes, dtype=int)
        self._lake_map.fill(self._grid.BAD_INDEX)

    @staticmethod
    def _validate_routing(routing):
        if routing not in {"D4", "D8"}:
            raise ValueError(
                f"routing method not understood ({routing} not one of"
                f" {', '.join(repr(s) for s in ('D4', 'D8'))})."
            )
        return routing

    @property
    def routing(self):
        """Routing method used when finding depressions."""
        return self._routing

    def updated_boundary_conditions(self):
        """Update neighbor nodes for changed boundary conditions."""
        if self.routing == "D8":
            neighbors = self.grid.d8_adjacent_nodes_at_node.copy()
        else:
            neighbors = self.grid.adjacent_nodes_at_node.copy()
        neighbors[self.grid.status_at_node[neighbors] == NodeStatus.CLOSED] = -1
        self._neighbor_nodes_at_node = neighbors

    @property
    def flood_status(self):
        """Map of flood status (PIT, CURRENT_LAKE, UNFLOODED, or FLOODED)."""
        return self._flood_status

    @property
    def receivers(self):
        """Return array indicating which node receives flow."""
        return self._receivers

    @receivers.setter
    def receivers(self, receivers):
        self._receivers = receivers

    @property
    def depression_depth(self):
        """Return array of depression depths."""
        return self._depression_depth

    @property
    def depression_outlet_map(self):
        """Return array indicating the node-id of the depression outlet."""
        return self._depression_outlet_map

    @staticmethod
    def find_pit_nodes(grid, value_at_node, routing=None):
        """Locate local depressions in a gridded elevation field.

        Parameters
        ----------
        grid : ModelGrid
            A Landlab grid.
        value_at_node : array_like
            Array of values overwhich to look for pits.
        routing : {'D4', 'D8'}, optional
            The routing method to use. For raster grids, 'D8' will
            consider all of a nodes neighbors (including the diagonals).
            Otherwise, only neighboring nodes connected by links are
            considered.

        Returns
        -------
        ndarray of bool
            A boolean array that indicates if a value in the input array
            is a pit.

        Notes
        -----

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
        value_at_node = np.asarray(value_at_node).reshape(-1)
        if len(value_at_node) != grid.number_of_nodes:
            raise ValueError(
                f"array length mismatch (expected {grid.number_of_nodes}"
                f" but got {len(value_at_node)})"
            )

        if routing is not None:
            routing = DepressionFinderAndRouter._validate_routing(routing)

        if isinstance(grid, RasterModelGrid):
            routing = routing if routing else "D8"
        else:
            warnings.warn(
                f"ignoring supplied routing method {routing!r}, not a raster grid",
                stacklevel=2,
            )

        is_pit = np.full(grid.number_of_nodes, False)
        is_pit[grid.status_at_node == NodeStatus.CORE] = True

        is_open_node = np.full(grid.number_of_nodes, False)
        is_open_node[grid.status_at_node == NodeStatus.FIXED_VALUE] = True

        if routing == "D8":
            links = grid.active_d8
            nodes_at_link = grid.nodes_at_d8
        else:
            links = grid.active_links
            nodes_at_link = grid.nodes_at_link

        find_pits(value_at_node, nodes_at_link, links, is_open_node, is_pit)

        return is_pit

    @staticmethod
    def _user_pits_to_array(grid, pits):
        """Convert pits to an n-node long boolean array.

        Parameters
        ----------
        grid : ModelGrid
            A landlab grid.
        pits : str or array_like
            If a string, the name of an *at-nodes* field. Otherwise, *pits*
            is interpreted as either a boolean array of nodes of the pits, or an
            array of pit node IDs.

        Returns
        -------
        ndarray of bool
            An n-nodes long array of booleans identifying pit nodes.

        Examples
        --------
        >>> from landlab import RasterModelGrid
        >>> from landlab.components.depression_finder import DepressionFinderAndRouter
        >>> grid = RasterModelGrid((3, 4))
        >>> pits = [
        ...     [False, False, False, False],
        ...     [False, True, False, False],
        ...     [False, False, False, False],
        ... ]

        >>> DepressionFinderAndRouter._user_pits_to_array(grid, pits).reshape(
        ...     grid.shape
        ... )
        array([[False, False, False, False],
               [False,  True, False, False],
               [False, False, False, False]])

        >>> grid.at_node["pits"] = pits
        >>> DepressionFinderAndRouter._user_pits_to_array(grid, "pits").reshape(
        ...     grid.shape
        ... )
        array([[False, False, False, False],
               [False,  True, False, False],
               [False, False, False, False]])

        >>> DepressionFinderAndRouter._user_pits_to_array(grid, [5]).reshape(grid.shape)
        array([[False, False, False, False],
               [False,  True, False, False],
               [False, False, False, False]])
        """
        assert pits is not None

        if isinstance(pits, str):
            try:
                is_pit = grid.at_node[pits].copy()
            except FieldError as err:
                raise ValueError(
                    f"grid is missing the provided at-node field ({pits!r} not one of"
                    f" {', '.join(repr(s) for s in grid.at_node)})."
                ) from err
        else:
            pits = np.asarray(pits).reshape(-1)
            if pits.size != grid.number_of_nodes:
                is_pit = np.full(grid.number_of_nodes, False)
                is_pit[pits] = True
            else:
                is_pit = pits.copy()

        is_pit[grid.status_at_node != NodeStatus.CORE] = False

        return is_pit

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
        if self.routing == "D8":
            diag_nbrs = self._grid.diagonal_adjacent_nodes_at_node[the_node]
        else:
            diag_nbrs = None

        return links, nbrs, diag_nbrs

    def assign_outlet_receiver(self, outlet_node):
        """Find drainage direction for an outlet that does not flow into its own lake.

        Parameters
        ----------
        outlet_node : int
            An outlet node.

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
                self._flood_status[nbr] != FloodStatus.CURRENT_LAKE
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
        if self.routing == "D8":
            for nbr in diag_nbrs:
                # Again, to pass this first hurdle, the neighbor must:
                #   * not be part of the current lake
                #   * have a surface (if flooded, WATER surface)
                #     lower than our outlet node;
                #   * not be a closed boundary
                if (
                    self._flood_status[nbr] != FloodStatus.CURRENT_LAKE
                    and (
                        (self._elev[nbr] + self._depression_depth[nbr])
                        < self._elev[receiver]
                    )
                    and self._grid.status_at_node[nbr] != self._grid.BC_NODE_IS_CLOSED
                ):
                    # Next test: is it the steepest downhill grad so far?
                    # If so, we've found a candidate.
                    grad = (node_elev - self._elev[nbr]) / np.sqrt(
                        self.grid.dx**2 + self.grid.dy**2
                    )
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
        nbrs = self._neighbor_nodes_at_node[the_node]
        not_bad = nbrs != self._grid.BAD_INDEX
        not_too_high = self._elev[nbrs] < self._elev[the_node]
        not_current_lake = np.not_equal(
            self._flood_status[nbrs], FloodStatus.CURRENT_LAKE
        )
        not_flooded = np.not_equal(self._flood_status[nbrs], FloodStatus.FLOODED)

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
            self._flood_status[n] = FloodStatus.FLOODED
            self._depression_depth[n] = self._elev[outlet_id] - self._elev[n]
            self._depression_outlet_map[n] = outlet_id
            self._lake_map[n] = pit_node
            self._pits_flooded += 1
            pit_node_where = np.searchsorted(self._pit_nodes, pit_node)
            self._unique_pits[pit_node_where] = True
        elif np.any(fresh_nodes):  # lake is bigger than one or more existing
            self._flood_status[n] = FloodStatus.FLOODED
            depth_this_lake = self._elev[outlet_id] - self._elev[n]
            self._depression_depth[n] = depth_this_lake
            self._depression_outlet_map[n] = outlet_id
            # ^these two will just get stamped over as needed
            subsumed_lakes = np.unique(self._lake_map[n])  # IDed by pit_node
            # the final entry is self._grid.BAD_INDEX
            subs_lakes_where = np.searchsorted(self._pit_nodes, subsumed_lakes[1:])
            pit_node_where = np.searchsorted(self._pit_nodes, pit_node)
            self._unique_pits[subs_lakes_where] = False
            self._unique_pits[pit_node_where] = True
            self._pits_flooded -= subsumed_lakes.size - 2
            # -1 for the self._grid.BAD_INDEX that must be present; another -1
            # because a single lake is just replaced by a new lake.
            self._lake_map[n] = pit_node
        else:  # lake is subsumed within an existing lake
            print(" eaten lake")
            assert np.all(np.equal(self._flood_status[n], FloodStatus.CURRENT_LAKE))
            self._flood_status[n] = FloodStatus.FLOODED

    def find_depression_from_pit(self, pit_node, reroute_flow=True):
        """Find the extent of the nodes that form a pit.

        Identify extent of depression/lake whose lowest point is the node
        pit_node (which is a itself a pit, a.k.a., closed depression).

        Parameters
        ----------
        pit_node : int
            The node that is the lowest point of a pit.
        """
        # Flag the pit as being CURRENT_LAKE (it's the first node in the
        # current lake)
        self._flood_status[pit_node] = FloodStatus.CURRENT_LAKE

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
                self._neighbor_nodes_at_node,
                self.flood_status,
                self._elev,
                nodes_this_depression,
                pit_count,
            )

            if lowest_node_on_perimeter == nodes_this_depression[0]:
                raise NoOutletError(pit_node, pit_count)

            # note this can return the supplied node, if - somehow - the
            # surrounding nodes are all self._grid.BAD_INDEX
            # I BELIEVE THE IS_VALID_OUTLET FN SHOULD ASSIGN FLOW DIR
            found_outlet = self.is_valid_outlet(lowest_node_on_perimeter)

            # If we haven't found an outlet, add lowest_node to the lake list
            # and flag it as being part of the current lake/depression
            if not found_outlet:
                nodes_this_depression[pit_count] = lowest_node_on_perimeter
                self._flood_status[lowest_node_on_perimeter] = FloodStatus.CURRENT_LAKE
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
            if count >= max_count:
                raise RuntimeError(
                    "DepressionFinderAndRouter reached maximum number of iterations"
                    " when finding depressions"
                )
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
        self._unique_pits = np.zeros_like(self._pit_nodes, dtype=bool)
        # debug_count = 0
        for pit_node in self._pit_nodes:
            if self._flood_status[pit_node] != FloodStatus.PIT:
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

        >>> df = DepressionFinderAndRouter(rg, reroute_flow=False, pits=None)
        >>> df.map_depressions()
        >>> print(df.display_depression_map())
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

        if self._pits is None:
            is_pit = DepressionFinderAndRouter.find_pit_nodes(
                self.grid, self._elev, routing=self.routing
            )
        else:
            is_pit = DepressionFinderAndRouter._user_pits_to_array(
                self.grid, self._pits
            )

        self._pit_nodes = as_id_array(np.nonzero(is_pit)[0])

        # Set up "lake code" array
        self._flood_status.fill(FloodStatus.UNFLOODED)
        self._flood_status[self._pit_nodes] = FloodStatus.PIT

        self._identify_depressions_and_outlets(self._reroute_flow)

        if self._reroute_flow and ("flow__receiver_node" in self._grid.at_node):
            self._receivers = self._grid.at_node["flow__receiver_node"]
            self._sinks = self._grid.at_node["flow__sink_flag"]
            self._grads = self._grid.at_node["topographic__steepest_slope"]
            self._links = self._grid.at_node["flow__link_to_receiver_node"]
            self._route_flow()
            self._reaccumulate_flow()

    def _find_unresolved_neighbors(self, nbrs, nbr_links, receivers):
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
        >>> df._find_unresolved_neighbors(nbrs, nbr_links, rcvr)
        (array([30, 21]), array([43, 35]))
        >>> nbrs = rg.diagonal_adjacent_nodes_at_node[22]
        >>> nbr_links = rg.d8s_at_node[22, 4:]
        >>> df._find_unresolved_neighbors(nbrs, nbr_links, rcvr)
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

        if not isinstance(self.grid, RasterModelGrid) or self.routing == "D4":
            adjacent_nodes_at_node = self.grid.adjacent_nodes_at_node
            links_at_node = self.grid.links_at_node
            length_of_link = self.grid.length_of_link
            n_regular_neighbors = adjacent_nodes_at_node.shape[1]
        else:
            adjacent_nodes_at_node = self.grid.d8_adjacent_nodes_at_node
            links_at_node = self.grid.d8s_at_node
            length_of_link = self.grid.length_of_d8
            n_regular_neighbors = 4

        # We must now iterate until we've taken care of all the nodes in the
        # lake. In each iteration, we:
        #  1 - find the unresolved neighbors of nodes being processed
        #  2 - point them toward the nodes being processed
        #  3 - place them on the nodes_to_proc_next list
        # We stop when there are no more nodes to process.
        #    Note that the nested looping will be slow, but could be sped up
        # by translating to cython.
        counter = 0  # counts # of times thru loop as fail-safe
        while nodes_being_processed:
            # Get unresolved "regular" neighbors of the current nodes
            for node in nodes_being_processed:
                # Get active and unresolved neighbors of node
                (unresolved_nodes, unresolved_links) = self._find_unresolved_neighbors(
                    adjacent_nodes_at_node[node, :n_regular_neighbors],
                    links_at_node[node, :n_regular_neighbors],
                    self._receivers,
                )

                # They will now flow to node
                if len(unresolved_nodes):
                    self._receivers[unresolved_nodes] = node
                    if "flow__link_to_receiver_node" in self._grid.at_node:
                        self._links[unresolved_nodes] = unresolved_links
                        slopes = (
                            self._elev[unresolved_nodes] - self._elev[node]
                        ) / length_of_link[unresolved_links]
                        self._grads[unresolved_nodes] = np.maximum(slopes, 0.0)

                # Place them on the list of nodes to process next
                nodes_to_proc_next.extend(unresolved_nodes)

            # If we're working with a raster that has diagonals, do the same
            # for the diagonal neighbors
            if self.routing == "D8":
                # Get unresolved "regular" neighbors of the current nodes
                for node in nodes_being_processed:
                    (unresolved_nodes, unresolved_links) = (
                        self._find_unresolved_neighbors(
                            adjacent_nodes_at_node[node, n_regular_neighbors:],
                            links_at_node[node, n_regular_neighbors:],
                            self._receivers,
                        )
                    )
                    # They will now flow to node
                    if len(unresolved_nodes):
                        self._receivers[unresolved_nodes] = node
                        if "flow__link_to_receiver_node" in self._grid.at_node:
                            self._links[unresolved_nodes] = unresolved_links
                            slopes = (
                                self._elev[unresolved_nodes] - self._elev[node]
                            ) / length_of_link[unresolved_links]
                            self._grads[unresolved_nodes] = np.maximum(slopes, 0.0)

                    # Place them on the list of nodes to process next
                    nodes_to_proc_next.extend(unresolved_nodes)

            # Move to the next set of nodes
            nodes_being_processed = nodes_to_proc_next
            nodes_to_proc_next = []

            # Just in case
            counter += 1
            if counter >= self._grid.number_of_nodes:
                raise RuntimeError(
                    "DepressionFinderAndRouter reached maximum number of iterations"
                    " when routing flow"
                )

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
                if self.routing == "D8":
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

        self._sinks[self._pit_nodes] = False

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

    def display_depression_map(self):
        """Print a simple character-based map of depressions/lakes."""
        # Find the outlet nodes (just for display purposes)
        is_outlet = np.full(self.grid.number_of_nodes, False)
        is_flooded = (self._flood_status == FloodStatus.FLOODED) & (
            self.grid.status_at_node == NodeStatus.CORE
        )
        is_outlet[self._depression_outlet_map[is_flooded]] = True

        def get_symbol(n):
            if is_outlet[n]:
                return "o"
            elif self._flood_status[n] == FloodStatus.UNFLOODED:
                return "."
            else:
                return "~"

        symbol_array = np.vectorize(get_symbol)(
            range(self.grid.number_of_nodes)
        ).reshape(self.grid.shape)

        stream = StringIO()
        np.savetxt(stream, symbol_array, fmt="%c")
        return stream.getvalue()

    @property
    def lake_outlets(self):
        """Return the *unique* outlets for each lake.

        Outlets are returned in same order as given from `lake_codes`.
        """
        return np.array(self._depression_outlets)[self._unique_pits]

    @property
    def lake_codes(self):
        """Return the *unique* code assigned to each unique lake.

        These are the values used to map the lakes in the property
        `lake_map`.
        """
        return self._pit_nodes[self._unique_pits]

    @property
    def number_of_lakes(self):
        """Return the number of individual lakes."""
        return self._unique_pits.sum()

    @property
    def lake_map(self):
        """Return unique codes for each lake.

        Return an array of ints, where each node within a lake is labelled
        with a unique (non-consecutive) code corresponding to each unique lake.

        The codes used can be obtained with `lake_codes`. Nodes not in a
        lake are labelled with `ModelGrid.BAD_INDEX`.
        """
        return self._lake_map

    @property
    def lake_at_node(self):
        """Return an array that indicates if a node is flooded.

        Return a ``bool`` array, ``True`` if the node is flooded, ``False``
        otherwise.
        """
        return self._lake_map != self._grid.BAD_INDEX

    @property
    def lake_areas(self):
        """Return an array of the area of each lake.

        The order is the same as that returned by `lake_codes`.
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
        """Return an array of the volume of each lake.

        The order is the same as that returned by `lake_codes`.
        """
        lake_vols = np.empty(self.number_of_lakes)
        col_vols = self._grid.cell_area_at_node * self._depression_depth
        for lake_counter, lake_code in enumerate(self.lake_codes):
            each_cell_in_lake = col_vols[self._lake_map == lake_code]
            lake_vols[lake_counter] = each_cell_in_lake.sum()
        return lake_vols

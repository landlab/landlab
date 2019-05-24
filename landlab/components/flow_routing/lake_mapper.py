# -*- coding: utf-8 -*-
"""Find depressions on a topographic surface.

.. codeauthor:: gtucker, DEJH (Flow routing)
"""
# Routing by DEJH, Oct 15.
from __future__ import print_function

import numpy as np

import landlab
from landlab import (
    CLOSED_BOUNDARY,
    CORE_NODE,
    FIXED_VALUE_BOUNDARY,
    Component,
    FieldError,
    ModelParameterDictionary,
    RasterModelGrid,
)
from landlab.components.flow_accum import flow_accum_bw
from landlab.core.messages import error_message, warning_message
from landlab.core.model_parameter_dictionary import MissingKeyError
from landlab.core.utils import as_id_array
from landlab.grid.base import BAD_INDEX_VALUE as LOCAL_BAD_INDEX_VALUE

# Codes for depression status
_UNFLOODED = 0
_PIT = 1
_CURRENT_LAKE = 2
_FLOODED = 3

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

    Note the routing part of this component is not yet compatible with
    irregular grids.

    The prinary method of this class is
    *map_depressions(pits='flow__sink_flag', reroute_flow=True)*.

    Examples
    --------
    Route flow across a depression in a sloped surface.

    >>> from landlab import RasterModelGrid
    >>> from landlab.components import FlowAccumulator, DepressionFinderAndRouter
    >>> mg = RasterModelGrid((7, 7), xy_spacing=0.5)
    >>> z = mg.add_field('node', 'topographic__elevation', mg.node_x.copy())
    >>> z += 0.01 * mg.node_y
    >>> mg.at_node['topographic__elevation'].reshape(mg.shape)[2:5, 2:5] *= 0.1
    >>> fr = FlowAccumulator(mg, flow_director='D8')
    >>> fr.run_one_step()  # the flow "gets stuck" in the hole
    >>> mg.at_node['flow__receiver_node'].reshape(mg.shape)
    array([[ 0,  1,  2,  3,  4,  5,  6],
           [ 7,  7, 16, 17, 18, 18, 13],
           [14, 14, 16, 16, 17, 18, 20],
           [21, 21, 16, 23, 24, 25, 27],
           [28, 28, 23, 30, 31, 32, 34],
           [35, 35, 30, 31, 32, 32, 41],
           [42, 43, 44, 45, 46, 47, 48]])
    >>> mg.at_node['drainage_area'].reshape(mg.shape)
    array([[ 0.  ,  0.  ,  0.  ,  0.  ,  0.  ,  0.  ,  0.  ],
           [ 0.25,  0.25,  0.25,  0.25,  0.25,  0.25,  0.  ],
           [ 0.25,  0.25,  5.  ,  1.5 ,  1.  ,  0.25,  0.  ],
           [ 0.25,  0.25,  3.  ,  0.75,  0.5 ,  0.25,  0.  ],
           [ 0.25,  0.25,  2.  ,  1.5 ,  1.  ,  0.25,  0.  ],
           [ 0.25,  0.25,  0.25,  0.25,  0.25,  0.25,  0.  ],
           [ 0.  ,  0.  ,  0.  ,  0.  ,  0.  ,  0.  ,  0.  ]])
    >>> df = DepressionFinderAndRouter(mg)
    >>> df.map_depressions()  # reroute_flow defaults to True
    >>> mg.at_node['flow__receiver_node'].reshape(mg.shape)
    array([[ 0,  1,  2,  3,  4,  5,  6],
           [ 7,  7, 16, 17, 18, 18, 13],
           [14, 14,  8, 16, 17, 18, 20],
           [21, 21, 16, 16, 24, 25, 27],
           [28, 28, 23, 24, 24, 32, 34],
           [35, 35, 30, 31, 32, 32, 41],
           [42, 43, 44, 45, 46, 47, 48]])
    >>> mg.at_node['drainage_area'].reshape(mg.shape)
    array([[ 0.  ,  0.  ,  0.  ,  0.  ,  0.  ,  0.  ,  0.  ],
           [ 5.25,  5.25,  0.25,  0.25,  0.25,  0.25,  0.  ],
           [ 0.25,  0.25,  5.  ,  1.5 ,  1.  ,  0.25,  0.  ],
           [ 0.25,  0.25,  0.75,  2.25,  0.5 ,  0.25,  0.  ],
           [ 0.25,  0.25,  0.5 ,  0.5 ,  1.  ,  0.25,  0.  ],
           [ 0.25,  0.25,  0.25,  0.25,  0.25,  0.25,  0.  ],
           [ 0.  ,  0.  ,  0.  ,  0.  ,  0.  ,  0.  ,  0.  ]])
    >>> df.lake_at_node.reshape(mg.shape)  # doctest: +NORMALIZE_WHITESPACE
    array([[False, False, False, False, False, False, False],
           [False, False, False, False, False, False, False],
           [False, False,  True,  True,  True, False, False],
           [False, False,  True,  True,  True, False, False],
           [False, False,  True,  True,  True, False, False],
           [False, False, False, False, False, False, False],
           [False, False, False, False, False, False, False]], dtype=bool)
    >>> df.lake_map.reshape(mg.shape)  # doctest: +NORMALIZE_WHITESPACE
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
    array([ 2.25])

    Because rereoute_flow defaults to True, the flow connectivity fields
    created by the FlowAccumulator will have now been modified to route flow over
    the depressions in the surface. The topogrphy itself is not modified.
    """

    _name = "DepressionFinderAndRouter"

    _input_var_names = ("topographic__elevation",)

    _output_var_names = ("depression__depth", "depression__outlet_node")

    _var_units = {
        "topographic__elevation": "m",
        "depression__depth": "m",
        "depression__outlet_node": "-",
    }

    _var_mapping = {
        "topographic__elevation": "node",
        "depression__depth": "node",
        "depression__outlet_node": "node",
    }

    _var_doc = {
        "topographic__elevation": "Surface topographic elevation",
        "depression__depth": "Depth of depression below its spillway point",
        "depression__outlet_node": "If a depression, the id of the outlet node for that depression, "
        "otherwise BAD_INDEX_VALUE",
    }

    def __init__(self, grid, routing="D8"):
        """Create a DepressionFinderAndRouter.

        Constructor assigns a copy of the grid, sets the current time, and
        calls the initialize method.

        Parameters
        ----------
        grid : RasterModelGrid
            A landlab RasterModelGrid.
        routing : 'D8' or 'D4' (optional)
            If grid is a raster type, controls whether lake connectivity can
            occur on diagonals ('D8', default), or only orthogonally ('D4').
            Has no effect if grid is not a raster.
        """
        super(DepressionFinderAndRouter, self).__init__(grid)
        self._grid = grid
        self._bc_set_code = self.grid.bc_set_code

        if routing != "D8":
            assert routing == "D4"
        self._routing = routing

        if isinstance(grid, RasterModelGrid) and (routing == "D8"):
            self._D8 = True
            self.num_nbrs = 8
            self._diag_link_length = np.sqrt(grid._dx ** 2 + grid._dy ** 2)
        else:
            self._D8 = False  # useful shorthand for thia test we do a lot
            if isinstance(grid, RasterModelGrid):
                self.num_nbrs = 4
            else:
                self.num_nbrs = self.grid.links_at_node.shape[1]

        if "flow__receiver_node" in self._grid.at_node:
            if self._grid.at_node["flow__receiver_node"].size != self._grid.size(
                "node"
            ):
                msg = (
                    "A route-to-multiple flow director has been "
                    "run on this grid. The depression finder is "
                    "not compatible with the grid anymore. Use "
                    "DepressionFinderAndRouter with reroute_flow="
                    "True only with route-to-one methods. If using this "
                    "component with such a flow directing method is desired "
                    "please open a GitHub Issue/"
                )
                raise NotImplementedError(msg)

        self._initialize()

    def _initialize(self, input_stream=None):
        """Initialize the component from an input file.

        The BMI-style initialize method takes an optional input_stream
        parameter, which may be either a ModelParameterDictionary object or
        an input stream from which a ModelParameterDictionary can read values.

        Parameters
        ----------
        input_stream : str, file_like, or ModelParameterDictionary, optional
            ModelParameterDictionary that holds the input parameters.
        """
        # Create a ModelParameterDictionary for the inputs
        if input_stream is None:
            inputs = None
        elif isinstance(input_stream, ModelParameterDictionary):
            inputs = input_stream
        else:
            inputs = ModelParameterDictionary(input_stream)

        # Make sure the grid includes elevation data. This means either:
        #  1. The grid has a node field called 'topographic__elevation', or
        #  2. The input file has an item called 'ELEVATION_FIELD_NAME' *and*
        #     a field by this name exists in the grid.
        try:
            self._elev = self._grid.at_node["topographic__elevation"]
        except FieldError:
            try:
                topo_field_name = inputs.read_string("ELEVATION_FIELD_NAME")
            except AttributeError:
                error_message(
                    """Because your grid does not have a node field
                    called "topographic__elevation", you need to pass the
                    name of a text input file or ModelParameterDictionary,
                    and this file or dictionary needs to include the name
                    of another field in your grid that contains your
                    elevation data."""
                )
                raise
            except MissingKeyError:
                error_message(
                    """Because your grid does not have a node field
                    called "topographic__elevation", your input file (or
                    ModelParameterDictionary) must include an entry with
                    the key "ELEVATION_FIELD_NAME", which gives the name
                    of a field in your grid that contains your elevation
                    data."""
                )
                raise
            try:
                self._elev = self._grid.at_node[topo_field_name]
            except AttributeError:
                warning_message(
                    """Your grid does not seem to have a node field
                    called {0}""".format(
                        topo_field_name
                    )
                )

        # Create output variables.
        #
        # Note that we initialize depression
        # outlet ID to LOCAL_BAD_INDEX_VALUE (which is a major clue!)
        self.depression_depth = self._grid.add_zeros(
            "node", "depression__depth", noclobber=False
        )
        self.depression_outlet_map = self._grid.add_zeros(
            "node", "depression__outlet_node", dtype=int, noclobber=False
        )
        self.depression_outlet_map += LOCAL_BAD_INDEX_VALUE

        # Later on, we'll need a number that's guaranteed to be larger than the
        # highest elevation in the grid.
        self._BIG_ELEV = 1.0e99

        self.updated_boundary_conditions()

        self._lake_outlets = []  # a list of each unique lake outlet
        # ^note this is nlakes-long

        self.is_pit = self._grid.add_ones("node", "is_pit", dtype=bool, noclobber=False)
        self.flood_status = self._grid.add_zeros(
            "node", "flood_status_code", dtype=int, noclobber=False
        )
        self._lake_map = np.empty(self._grid.number_of_nodes, dtype=int)
        self._lake_map.fill(LOCAL_BAD_INDEX_VALUE)

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
            diag_nbrs[self.grid.status_at_node[diag_nbrs] == CLOSED_BOUNDARY] = -1
            self._node_nbrs = np.concatenate((self._node_nbrs, diag_nbrs), 1)
            self._link_lengths = np.empty(8, dtype=float)
            self._link_lengths[0] = dx
            self._link_lengths[2] = dx
            self._link_lengths[1] = dy
            self._link_lengths[3] = dy
            self._link_lengths[4:].fill(np.sqrt(dx * dx + dy * dy))
        elif (type(self.grid) is landlab.grid.raster.RasterModelGrid) and (
            self._routing == "D4"
        ):
            self._link_lengths = np.empty(4, dtype=float)
            self._link_lengths[0] = dx
            self._link_lengths[2] = dx
            self._link_lengths[1] = dy
            self._link_lengths[3] = dy
        else:
            self._link_lengths = self.grid.length_of_link

    def _find_pits(self):
        """Locate local depressions ("pits") in a gridded elevation field.

        Notes
        -----
        **Uses**:

        * ``self._elev``
        * ``self._grid``

        **Creates**:

        * ``self.is_pit`` (node array of booleans): Flag indicating whether
          the node is a pit.
        * ``self.number_of_pits`` (int): Number of pits found.
        * ``self.pit_node_ids`` (node array of ints): IDs of the nodes that
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
        self.is_pit.fill(True)
        self.is_pit[self._grid.boundary_nodes] = False

        # Loop over all active links: if one of a link's two nodes is higher
        # than the other, the higher one is not a pit. Also, if they have
        # equal elevations and one is an open boundary, the other is not a pit.
        act_links = self._grid.active_links
        h_orth = self._grid.node_at_link_head[act_links]
        t_orth = self._grid.node_at_link_tail[act_links]

        # These two lines assign the False flag to any node that is higher
        # than its partner on the other end of its link
        self.is_pit[
            h_orth[np.where(self._elev[h_orth] > self._elev[t_orth])[0]]
        ] = False
        self.is_pit[
            t_orth[np.where(self._elev[t_orth] > self._elev[h_orth])[0]]
        ] = False

        # If we have a raster grid, handle the diagonal active links too
        # (At the moment, their data structure is a bit different)
        # TODO: update the diagonal link data structures
        # DEJH doesn't understand why this can't be vectorized as above...
        if self._D8:
            for h, t in self.grid.nodes_at_diagonal[self.grid.active_diagonals]:
                if self._elev[h] > self._elev[t]:
                    self.is_pit[h] = False
                elif self._elev[t] > self._elev[h]:
                    self.is_pit[t] = False
                elif self._elev[h] == self._elev[t]:
                    if self._grid.status_at_node[h] == FIXED_VALUE_BOUNDARY:
                        self.is_pit[t] = False
                    elif self._grid.status_at_node[t] == FIXED_VALUE_BOUNDARY:
                        self.is_pit[h] = False

        # Record the number of pits and the IDs of pit nodes.
        self.number_of_pits = np.count_nonzero(self.is_pit)
        self.pit_node_ids = as_id_array(np.where(self.is_pit)[0])

    def find_lowest_node_on_lake_perimeter(self, nodes_this_depression):
        """Locate the lowest node on the margin of the "lake".

        Parameters
        ----------
        nodes_this_depression : array_like of int
            Nodes that form a pit.

        Returns
        -------
        int
            The lowest node on the perimeter of a depression.
        """
        # Start with the first node on the list, and an arbitrarily large elev
        lowest_node = nodes_this_depression[0]
        lowest_elev = self._BIG_ELEV

        for n in nodes_this_depression:

            for nbr in self._node_nbrs[n]:
                if nbr != -1:
                    if self.flood_status[nbr] == _UNFLOODED:
                        if self._elev[nbr] < lowest_elev:
                            lowest_node = nbr
                            lowest_elev = self._elev[nbr]
                    elif (
                        self.flood_status[nbr] == _PIT
                        or self.flood_status[nbr] == _FLOODED
                    ):
                        nodes_this_depression.append(nbr)
                        self.flood_status[nbr] = _CURRENT_LAKE
        if lowest_elev == self._BIG_ELEV:
            print("Unable to find drainage outlet for a lake.")
            print("In lake with " + str(len(nodes_this_depression)), "nodes:")
            print(str(nodes_this_depression))

            for i in nodes_this_depression:
                print("Node ID: ", i)
                print("Node Elevation: ", self._elev[i])
                print("Node Flood Status: ", self.flood_status[i])
                print("Node Grid Status: ", self._grid.status_at_node[i])
                print("Node Neigbors: ", self._node_nbrs[i])
                print("Neighbor Elevations: ", self._elev[self._node_nbrs[i]])
                print("Neigbor Flood Status: ", self.flood_status[self._node_nbrs[i]])
                print("Neigbor Status: ", self._grid.status_at_node[self._node_nbrs[i]])
            warning_message(
                """If you see no data values in any of the elevation terms
                this may because you have disconnected open nodes (which
                sometimes occurs during raster clipping.

                Consider running
                set_open_nodes_disconnected_from_watershed_to_closed
                which will remove isolated open nodes."""
            )

        assert lowest_elev < self._BIG_ELEV, "failed to find lowest perim node"
        return lowest_node

    def _links_and_nbrs_at_node(self, the_node):
        """Compile and return arrays with IDs of neighbor links and nodes.

        If D8 Raster, returns *diag_nbrs* containing the diagonal neighbors;
        otherwise, *diag_nbrs = None*.

        Examples
        --------
        >>> from landlab import RasterModelGrid
        >>> from landlab.components import DepressionFinderAndRouter
        >>> rg = RasterModelGrid((3, 3))
        >>> z = rg.add_zeros('node', 'topographic__elevation')
        >>> z[4] = 2.0
        >>> df = DepressionFinderAndRouter(rg, routing='D4')
        >>> (links, nbrs, dnbrs) = df._links_and_nbrs_at_node(4)
        >>> links
        array([6, 8, 5, 3])
        >>> nbrs
        array([5, 7, 3, 1])
        >>> dnbrs
        >>> df = DepressionFinderAndRouter(rg, routing='D8')
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
        >>> rg.status_at_node[rg.nodes_at_right_edge] = CLOSED_BOUNDARY
        >>> z = rg.add_zeros('node', 'topographic__elevation')
        >>> z[:] = rg.x_of_node + 0.01 * rg.y_of_node
        >>> lake_nodes = np.array([10, 16, 17, 18, 24, 32, 33, 38, 40])
        >>> z[lake_nodes] *= 0.1
        >>> fr = FlowAccumulator(rg, flow_director='D4')
        >>> fr.run_one_step()
        >>> rg.at_node['flow__receiver_node']
        array([ 0,  1,  2,  3,  4,  5,  6,  7,  7, 16, 10, 10, 11, 13, 14, 14, 16,
               16, 17, 18, 20, 21, 21, 16, 17, 24, 33, 27, 28, 28, 29, 24, 32, 32,
               34, 35, 35, 38, 38, 38, 33, 41, 42, 43, 44, 45, 46, 47, 48])
        >>> df = DepressionFinderAndRouter(rg, routing='D4')
        >>> df.map_depressions()
        >>> rg.at_node['flow__receiver_node']
        array([ 0,  1,  2,  3,  4,  5,  6,  7,  7, 16, 17, 10, 11, 13, 14, 14, 15,
               16, 17, 18, 20, 21, 21, 16, 17, 24, 33, 27, 28, 28, 29, 38, 31, 32,
               34, 35, 35, 36, 37, 38, 33, 41, 42, 43, 44, 45, 46, 47, 48])
        >>> fr = FlowAccumulator(rg, flow_director='D8')
        >>> fr.run_one_step()
        >>> rg.at_node['flow__receiver_node']
        array([ 0,  1,  2,  3,  4,  5,  6,  7,  7, 16, 16, 10, 18, 13, 14, 14, 16,
               16, 17, 18, 20, 21, 21, 16, 16, 24, 33, 27, 28, 28, 24, 24, 24, 32,
               34, 35, 35, 38, 38, 38, 32, 41, 42, 43, 44, 45, 46, 47, 48])
        >>> df = DepressionFinderAndRouter(rg, routing='D8')
        >>> df.map_depressions()
        >>> rg.at_node['flow__receiver_node']
        array([ 0,  1,  2,  3,  4,  5,  6,  7,  7, 16, 16, 10, 18, 13, 14, 14,  8,
               16, 17, 18, 20, 21, 21, 16, 16, 24, 33, 27, 28, 28, 24, 24, 24, 32,
               34, 35, 35, 38, 32, 38, 32, 41, 42, 43, 44, 45, 46, 47, 48])
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
                self.flood_status[nbr] != _CURRENT_LAKE
                and (
                    (self._elev[nbr] + self.depression_depth[nbr])
                    < self._elev[receiver]
                )
                and self._grid.status_at_node[nbr] != CLOSED_BOUNDARY
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
                    self.flood_status[nbr] != _CURRENT_LAKE
                    and (
                        (self._elev[nbr] + self.depression_depth[nbr])
                        < self._elev[receiver]
                    )
                    and self._grid.status_at_node[nbr] != CLOSED_BOUNDARY
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
        not_bad = nbrs != LOCAL_BAD_INDEX_VALUE
        not_too_high = self._elev[nbrs] < self._elev[the_node]
        not_current_lake = np.not_equal(self.flood_status[nbrs], _CURRENT_LAKE)
        not_flooded = np.not_equal(self.flood_status[nbrs], _FLOODED)

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
                    dep_out = self.depression_outlet_map[nbrs[i]]
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
        if np.any(all_probs):
            return True
        else:
            return False

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
        if self._grid.status_at_node[the_node] == FIXED_VALUE_BOUNDARY:
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
        fresh_nodes = np.equal(self._lake_map[n], LOCAL_BAD_INDEX_VALUE)
        if np.all(fresh_nodes):  # a new lake
            self.flood_status[n] = _FLOODED
            self.depression_depth[n] = self._elev[outlet_id] - self._elev[n]
            self.depression_outlet_map[n] = outlet_id
            self._lake_map[n] = pit_node
            self._pits_flooded += 1
            pit_node_where = np.searchsorted(self.pit_node_ids, pit_node)
            self._unique_pits[pit_node_where] = True
        elif np.any(fresh_nodes):  # lake is bigger than one or more existing
            self.flood_status[n] = _FLOODED
            depth_this_lake = self._elev[outlet_id] - self._elev[n]
            self.depression_depth[n] = depth_this_lake
            self.depression_outlet_map[n] = outlet_id
            # ^these two will just get stamped over as needed
            subsumed_lakes = np.unique(self._lake_map[n])  # IDed by pit_node
            # the final entry is LOCAL_BAD_INDEX_VALUE
            subs_lakes_where = np.searchsorted(self.pit_node_ids, subsumed_lakes[1:])
            pit_node_where = np.searchsorted(self.pit_node_ids, pit_node)
            self._unique_pits[subs_lakes_where] = False
            self._unique_pits[pit_node_where] = True
            self._pits_flooded -= subsumed_lakes.size - 2
            # -1 for the LOCAL_BAD_INDEX_VALUE that must be present; another -1
            # because a single lake is just replaced by a new lake.
            self._lake_map[n] = pit_node
        else:  # lake is subsumed within an existing lake
            print(" eaten lake")
            assert np.all(np.equal(self.flood_status[n], _CURRENT_LAKE))
            self.flood_status[n] = _FLOODED

    def find_depression_from_pit(self, pit_node, reroute_flow=True):
        """Find the extent of the nodes that form a pit.

        Identify extent of depression/lake whose lowest point is the node
        pit_node (which is a itself a pit, a.k.a., closed depression).

        Parameters
        ----------
        pit_node : int
            The node that is the lowest point of a pit.
        """

        # Place pit_node at top of depression list
        nodes_this_depression = []
        nodes_this_depression.insert(0, pit_node)

        # Flag the pit as being _CURRENT_LAKE (it's the first node in the
        # current lake)
        self.flood_status[pit_node] = _CURRENT_LAKE

        # This flag keeps track of when we're done with this depression
        found_outlet = False

        # Safety check
        count = 0
        max_count = self._grid.number_of_nodes + 1

        while not found_outlet:
            lowest_node_on_perimeter = self.find_lowest_node_on_lake_perimeter(
                nodes_this_depression
            )
            # note this can return the supplied node, if - somehow - the
            # surrounding nodes are all LOCAL_BAD_INDEX_VALUE
            # I BELIEVE THE IS_VALID_OUTLET FN SHOULD ASSIGN FLOW DIR
            found_outlet = self.is_valid_outlet(lowest_node_on_perimeter)

            # If we haven't found an outlet, add lowest_node to the lake list
            # and flag it as being part of the current lake/depression
            if not found_outlet:
                nodes_this_depression.append(lowest_node_on_perimeter)
                self.flood_status[lowest_node_on_perimeter] = _CURRENT_LAKE

            # If we HAVE found an outlet, and we are re-routing flow, then
            # assign the proper flow direction to the outlet node. If it is an
            # open boundary, then it drains to itself. Otherwise, call
            # assign_outlet_receiver to find the correct receiver (so that it
            # doesn't simply drain back into the lake)
            elif ("flow__receiver_node" in self._grid.at_node) and reroute_flow:
                if self._grid.status_at_node[lowest_node_on_perimeter] != CORE_NODE:
                    self._grid.at_node["flow__receiver_node"][
                        lowest_node_on_perimeter
                    ] = lowest_node_on_perimeter
                else:
                    self.assign_outlet_receiver(lowest_node_on_perimeter)

            # Safety check, in case a bug (ha!) puts us in an infinite loop
            assert count < max_count, "too many iterations in lake filler!"
            count += 1

        self.depression_outlets.append(lowest_node_on_perimeter)
        # Now that we've mapped this depression, record it in the arrays
        # depression_depth, depression_outlet, and flood_status
        self._record_depression_depth_and_outlet(
            nodes_this_depression, lowest_node_on_perimeter, pit_node
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
        self._unique_pits = np.zeros_like(self.pit_node_ids, dtype=bool)
        # debug_count = 0
        for pit_node in self.pit_node_ids:
            if self.flood_status[pit_node] != _PIT:
                from landlab import BAD_INDEX_VALUE

                self.depression_outlets.append(BAD_INDEX_VALUE)
            else:
                self.find_depression_from_pit(pit_node, reroute_flow)
                self._pits_flooded += 1

        assert len(self.depression_outlets) == self._unique_pits.size

        self.unique_lake_outlets = np.array(self.depression_outlets)[self._unique_pits]

    def map_depressions(self, pits="flow__sink_flag", reroute_flow=True):
        """Map depressions/lakes in a topographic surface.

        Parameters
        ----------
        pits : array or str or None, optional
            If a field name, the boolean field containing True where pits.
            If an array, either a boolean array of nodes of the pits, or an
            array of pit node IDs. It does not matter whether or not open
            boundary nodes are flagged as pits; they are never treated as such.
            Default is 'flow__sink_flag', the pit field output from
            'route_flow_dn'
        reroute_flow : bool, optional
            If True (default), and the component detects the output fields in
            the grid produced by the route_flow_dn component, this component
            will modify the existing flow fields to route the flow across the
            lake surface(s).
            Ensure you call this method *after* you have already routed flow
            in each loop of your model.

        Examples
        --------
        Test #1: 5x5 raster grid with a diagonal lake.

        >>> import numpy as np
        >>> from landlab import RasterModelGrid
        >>> from landlab.components.flow_routing import (
        ...     DepressionFinderAndRouter)

        >>> rg = RasterModelGrid((5, 5))
        >>> z = rg.add_zeros('node', 'topographic__elevation')
        >>> z[:] = np.array([100., 100.,  95., 100., 100.,
        ...                  100., 101.,  92.,   1., 100.,
        ...                  100., 101.,   2., 101., 100.,
        ...                  100.,   3., 101., 101., 100.,
        ...                   90.,  95., 100., 100., 100.])
        >>> df = DepressionFinderAndRouter(rg)
        >>> df.map_depressions(pits=None, reroute_flow=False)
        >>> df.display_depression_map()  # doctest: +NORMALIZE_WHITESPACE
        . . . . .
        . . . ~ .
        . . ~ . .
        . ~ . . .
        o . . . .
        """
        if self._bc_set_code != self.grid.bc_set_code:
            self.updated_boundary_conditions()
            self._bc_set_code = self.grid.bc_set_code
        self._lake_map.fill(LOCAL_BAD_INDEX_VALUE)
        self.depression_outlet_map.fill(LOCAL_BAD_INDEX_VALUE)
        self.depression_depth.fill(0.0)
        self.depression_outlets = []  # reset these
        # Locate nodes with pits
        if type(pits) == str:
            try:
                pits = self._grid.at_node[pits]
                supplied_pits = np.where(pits)[0]
                self.pit_node_ids = as_id_array(
                    np.setdiff1d(supplied_pits, self._grid.boundary_nodes)
                )
                self.number_of_pits = self.pit_node_ids.size
                self.is_pit.fill(False)
                self.is_pit[self.pit_node_ids] = True
            except FieldError:
                self._find_pits()
        elif pits is None:
            self._find_pits()
        else:  # hopefully an array or other sensible iterable
            if len(pits) == self._grid.number_of_nodes:
                supplied_pits = np.where(pits)[0]
            else:  # it's an array of node ids
                supplied_pits = pits
            # remove any boundary nodes from the supplied pit list
            self.pit_node_ids = as_id_array(
                np.setdiff1d(supplied_pits, self._grid.boundary_nodes)
            )

            self.number_of_pits = self.pit_node_ids.size
            self.is_pit.fill(False)
            self.is_pit[self.pit_node_ids] = True
        # Set up "lake code" array
        self.flood_status.fill(_UNFLOODED)
        self.flood_status[self.pit_node_ids] = _PIT

        self._identify_depressions_and_outlets(reroute_flow)

        if reroute_flow and ("flow__receiver_node" in self._grid.at_node):

            self.receivers = self._grid.at_node["flow__receiver_node"]
            self.sinks = self._grid.at_node["flow__sink_flag"]
            self.grads = self._grid.at_node["topographic__steepest_slope"]
            self.links = self._grid.at_node["flow__link_to_receiver_node"]
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
        >>> z = rg.add_zeros('node', 'topographic__elevation')
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
        >>> z = rg.add_zeros('node', 'topographic__elevation')
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
        >>> z = rg.add_zeros('node', 'topographic__elevation')
        >>> rcvr = rg.add_zeros('node', 'flow__receiver_node', dtype=int)
        >>> rcvr[:] = np.arange(rg.number_of_nodes)
        >>> lake_nodes = np.array([10, 12, 13, 19, 20, 21, 25, 26, 27, 28, 29, 30, 33, 34, 35, 36, 37, 38, 44, 45, 46])
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
        self.receivers[lake_nodes] = UNRESOLVED

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
                    self.grid.adjacent_nodes_at_node[cn],
                    self.grid.links_at_node[cn],
                    self.receivers,
                )
                # They will now flow to cn
                if nbrs.size > 0:
                    self.receivers[nbrs] = cn
                    if "flow__link_to_receiver_node" in self._grid.at_node:
                        self.links[nbrs] = lnks
                        slopes = (
                            self._elev[nbrs] - self._elev[cn]
                        ) / self._grid.length_of_link[lnks]
                        self.grads[nbrs] = np.maximum(slopes, 0.0)

                # Place them on the list of nodes to process next
                for n in nbrs:
                    nodes_to_proc_next.append(n)

            # If we're working with a raster that has diagonals, do the same
            # for the diagonal neighbors
            if self._D8:

                # Get unresolved "regular" neighbors of the current nodes
                for cn in nodes_being_processed:

                    # Get active and unresolved diagonal neighbors of cn
                    #                    nbrs = self._find_unresolved_neighbors(
                    #                            self._grid._get_diagonal_list(cn), self.receivers)
                    (nbrs, diags) = self._find_unresolved_neighbors_new(
                        self._grid.diagonal_adjacent_nodes_at_node[cn],
                        self._grid.d8s_at_node[cn, 4:],
                        self.receivers,
                    )

                    # They will now flow to cn
                    if nbrs.size > 0:
                        self.receivers[nbrs] = cn
                        if "flow__link_to_receiver_node" in self._grid.at_node:
                            self.links[nbrs] = diags
                            slopes = (
                                self._elev[nbrs] - self._elev[cn]
                            ) / self._diag_link_length
                            self.grads[nbrs] = np.maximum(slopes, 0.0)

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
            nodes_in_lake = np.where(self.lake_map == lake_code)[0]
            if len(nodes_in_lake) > 0:

                # find the correct outlet for the lake, if necessary
                if self.lake_map[self.receivers[outlet_node]] == lake_code:
                    nbrs = self.grid.active_adjacent_nodes_at_node[outlet_node]
                    not_lake = nbrs[np.where(self.lake_map[nbrs] != lake_code)[0]]
                    min_index = np.argmin(self._elev[not_lake])
                    new_receiver = not_lake[min_index]

                    # set receiver for new outlet.
                    self.receivers[outlet_node] = new_receiver

                # reset_link for new outlet
                outlet_receiver = self.receivers[outlet_node]
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
                    new_link = LOCAL_BAD_INDEX_VALUE
                self.links[outlet_node] = new_link

                # make a check
                assert (
                    self.lake_map[self.receivers[outlet_node]] != lake_code
                ), "outlet of lake drains to itself!"

                # Route flow
                self._route_flow_for_one_lake(outlet_node, nodes_in_lake)

        self.sinks[self.pit_node_ids] = False

    def _reaccumulate_flow(self):
        """Update drainage area, discharge, upstream order, and flow link.

        Invoke the accumulator a second time to update drainage area,
        discharge, and upstream order.
        """
        # Calculate drainage area, discharge, and downstr->upstr order
        Q_in = self._grid.at_node["water__unit_flux_in"]
        areas = self._grid.cell_area_at_node.copy()
        areas[self._grid.closed_boundary_nodes] = 0.0

        self.a, q, s = flow_accum_bw.flow_accumulation(
            self.receivers, node_cell_area=areas, runoff_rate=Q_in
        )

        # finish the property updating:
        self.grid.at_node["drainage_area"][:] = self.a
        self.grid.at_node["surface_water__discharge"][:] = q
        self.grid.at_node["flow__upstream_node_order"][:] = s

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
            if isinstance(self._grid, landlab.grid.raster.RasterModelGrid):
                link_l = self._link_lengths
            else:  # Voronoi
                link_l = self._link_lengths[self._grid.links_at_node[outlet_node, :]]
            eff_slopes = (self._elev[outlet_node] - self._elev[out_draining]) / link_l[
                unique_indxs[1:]
            ]
            lowest = np.argmax(eff_slopes)
            lowest_node = out_draining[lowest]
            # route the flow
            self.receivers[outlet_node] = lowest_node
        else:
            self.receivers[outlet_node] = outlet_node

    def display_depression_map(self):
        """Print a simple character-based map of depressions/lakes."""
        # Find the outlet nodes (just for display purposes)
        is_outlet = np.zeros(self._grid.number_of_nodes, dtype=bool)
        for i in self._grid.core_nodes:
            if self.flood_status[i] == _FLOODED:
                is_outlet[self.depression_outlet_map[i]] = True

        n = 0
        for r in range(self._grid.number_of_node_rows):
            for c in range(self._grid.number_of_node_columns):
                if is_outlet[n]:
                    print("o", end=" ")
                elif self.flood_status[n] == _UNFLOODED:
                    print(".", end=" ")
                else:
                    print("~", end=" ")
                n += 1
            print()

    @property
    def lake_outlets(self):
        """Returns the *unique* outlets for each lake, in same order as the
        return from lake_codes."""
        return np.array(self.depression_outlets)[self._unique_pits]

    @property
    def lake_codes(self):
        """Returns the *unique* code assigned to each unique lake.

        These are the values used to map the lakes in the property
        "lake_map".
        """
        return self.pit_node_ids[self._unique_pits]

    @property
    def number_of_lakes(self):
        """Return the number of individual lakes."""
        return self._unique_pits.sum()

    @property
    def lake_map(self):
        """Return an array of ints, where each node within a lake is labelled
        with a unique (non-consecutive) code corresponding to each unique lake.

        The codes used can be obtained with *lake_codes*. Nodes not in a
        lake are labelled with LOCAL_BAD_INDEX_VALUE.
        """
        return self._lake_map

    @property
    def lake_at_node(self):
        """Return a boolean array, True if the node is flooded, False
        otherwise."""
        return self._lake_map != LOCAL_BAD_INDEX_VALUE

    @property
    def lake_areas(self):
        """A nlakes-long array of the area of each lake.

        The order is the same as that returned by *lake_codes*.
        """
        lake_areas = np.empty(self.number_of_lakes)
        lake_counter = 0
        for lake_code in self.lake_codes:
            each_cell_in_lake = self._grid.cell_area_at_node[self.lake_map == lake_code]
            lake_areas[lake_counter] = each_cell_in_lake.sum()
            lake_counter += 1
        return lake_areas

    @property
    def lake_volumes(self):
        """A nlakes-long array of the volume of each lake.

        The order is the same as that returned by *lake_codes*.
        """
        lake_vols = np.empty(self.number_of_lakes)
        lake_counter = 0
        col_vols = self._grid.cell_area_at_node * self.depression_depth
        for lake_code in self.lake_codes:
            each_cell_in_lake = col_vols[self.lake_map == lake_code]
            lake_vols[lake_counter] = each_cell_in_lake.sum()
            lake_counter += 1
        return lake_vols

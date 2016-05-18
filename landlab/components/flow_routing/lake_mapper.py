# -*- coding: utf-8 -*-
"""Find depressions on a topographic surface.

.. codeauthor:: gtucker, DEJH (Flow routing)
"""
# Routing by DEJH, Oct 15.
from __future__ import print_function

import numpy as np
from landlab import (ModelParameterDictionary, Component, FieldError,
                     FIXED_VALUE_BOUNDARY)
from landlab.core.utils import as_id_array
from landlab.core.model_parameter_dictionary import MissingKeyError
from landlab.components.flow_accum import flow_accum_bw
from landlab.grid.base import BAD_INDEX_VALUE as LOCAL_BAD_INDEX_VALUE
# LOCAL_BAD_INDEX_VALUE = np.iinfo(np.int32).max
import landlab


# Codes for depression status
_UNFLOODED = 0
_PIT = 1
_CURRENT_LAKE = 2
_FLOODED = 3


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

    Construction::

        DepressionFinderAndRouter(grid, grid, routing='D8')

    Parameters
    ----------
    grid : RasterModelGrid
        A landlab RasterModelGrid.
    routing : {'D8', 'D4'} (optional)
        If grid is a raster type, controls whether lake connectivity can
        occur on diagonals ('D8', default), or only orthogonally ('D4').
        Has no effect if grid is not a raster.

    Examples
    --------
    Route flow across a depression in a sloped surface.

    >>> from landlab import RasterModelGrid
    >>> from landlab.components import FlowRouter, DepressionFinderAndRouter
    >>> mg = RasterModelGrid((7, 7), 0.5)
    >>> _ = mg.add_field('node', 'topographic__elevation', mg.node_x.copy())
    >>> mg.at_node['topographic__elevation'].reshape(mg.shape)[2:5, 2:5] = 0.
    >>> fr = FlowRouter(mg)
    >>> fr.run_one_step()  # the flow "gets stuck" in the hole
    >>> mg.at_node['drainage_area'].reshape(mg.shape)
    array([[ 0.  ,  0.  ,  0.  ,  0.  ,  0.  ,  0.  ,  0.  ],
           [ 0.25,  0.25,  0.25,  0.25,  0.25,  0.25,  0.  ],
           [ 0.25,  0.25,  0.5 ,  0.5 ,  1.  ,  0.25,  0.  ],
           [ 0.25,  0.25,  0.25,  0.25,  0.5 ,  0.25,  0.  ],
           [ 0.25,  0.25,  0.5 ,  0.5 ,  1.  ,  0.25,  0.  ],
           [ 0.25,  0.25,  0.25,  0.25,  0.25,  0.25,  0.  ],
           [ 0.  ,  0.  ,  0.  ,  0.  ,  0.  ,  0.  ,  0.  ]])
    >>> df = DepressionFinderAndRouter(mg)
    >>> df.map_depressions()  # reroute_flow defaults to True
    >>> mg.at_node['drainage_area'].reshape(mg.shape)
    array([[ 0.  ,  0.  ,  0.  ,  0.  ,  0.  ,  0.  ,  0.  ],
           [ 0.25,  0.25,  0.25,  0.25,  0.25,  0.25,  0.  ],
           [ 5.25,  5.25,  3.75,  2.  ,  1.  ,  0.25,  0.  ],
           [ 0.25,  0.25,  1.25,  1.25,  0.5 ,  0.25,  0.  ],
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
    array([15])
    >>> df.lake_areas  # the area of each lake in lake_codes
    array([ 2.25])

    Because rereoute_flow defaults to True, the flow connectivity fields
    created by the FlowRouter will have now been modified to route flow over
    the depressions in the surface. The topogrphy itself is not modified.
    """

    _name = 'DepressionFinderAndRouter'

    _input_var_names = ('topographic__elevation',
                        )

    _output_var_names = ('depression__depth',
                         'depression__outlet_node',
                         )

    _var_units = {'topographic__elevation': 'm',
                  'depression__depth': 'm',
                  'depression__outlet_node': '-'
                  }

    _var_mapping = {'topographic__elevation': 'node',
                    'depression__depth': 'node',
                    'depression__outlet_node': 'node',
                    }

    _var_doc = {
        'topographic__elevation': 'Surface topographic elevation',
        'depression__depth': 'Depth of depression below its spillway point',
        'depression__outlet_node':
            'If a depression, the id of the outlet node for that depression, '
            'otherwise BAD_INDEX_VALUE'
    }

    def __init__(self, grid, routing='D8', **kwds):
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
        self._grid = grid
        if routing is not 'D8':
            assert routing is 'D4'
        self._routing = routing
        if ((type(self._grid) is landlab.grid.raster.RasterModelGrid) and
                (routing is 'D8')):
            self._D8 = True
            self.num_nbrs = 8
        else:
            self._D8 = False  # useful shorthand for thia test we do a lot
            if type(self._grid) is landlab.grid.raster.RasterModelGrid:
                self.num_nbrs = 4
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
        elif type(input_stream) == ModelParameterDictionary:
            inputs = input_stream
        else:
            inputs = ModelParameterDictionary(input_stream)

        # Make sure the grid includes elevation data. This means either:
        #  1. The grid has a node field called 'topographic__elevation', or
        #  2. The input file has an item called 'ELEVATION_FIELD_NAME' *and*
        #     a field by this name exists in the grid.
        try:
            self._elev = self._grid.at_node['topographic__elevation']
        except FieldError:
            try:
                topo_field_name = inputs.read_string('ELEVATION_FIELD_NAME')
            except AttributeError:
                print('Error: Because your grid does not have a node field')
                print('called "topographic__elevation", you need to pass the')
                print('name of a text input file or ModelParameterDictionary,')
                print('and this file or dictionary needs to include the name')
                print('of another field in your grid that contains your')
                print('elevation data.')
                raise AttributeError
            except MissingKeyError:
                print('Error: Because your grid does not have a node field')
                print('called "topographic__elevation", your input file (or')
                print('ModelParameterDictionary) must include an entry with')
                print('the key "ELEVATION_FIELD_NAME", which gives the name')
                print('of a field in your grid that contains your elevation')
                print('data.')
                raise MissingKeyError('ELEVATION_FIELD_NAME')
            try:
                self._elev = self._grid.at_node[topo_field_name]
            except AttributeError:
                print('Your grid does not seem to have a node field called',
                      topo_field_name)

        # Create output variables.
        #
        # Note that we initialize depression
        # outlet ID to LOCAL_BAD_INDEX_VALUE (which is a major clue!)
        self.depression_depth = self._grid.add_zeros('node',
                                                     'depression__depth',
                                                     noclobber=False)
        self.depression_outlet_map = self._grid.add_zeros(
            'node', 'depression__outlet_node', dtype=int, noclobber=False)
        self.depression_outlet_map += LOCAL_BAD_INDEX_VALUE

        # Later on, we'll need a number that's guaranteed to be larger than the
        # highest elevation in the grid.
        self._BIG_ELEV = np.amax(self._elev) + 1

        # We'll also need a handy copy of the node neighbor lists
        # TODO: presently, this grid method seems to only exist for Raster
        # grids. We need it for *all* grids!
        self._node_nbrs = self._grid.active_neighbors_at_node()
        dx = self._grid.dx
        dy = self._grid.dy
        if self._D8:
            diag_nbrs = self._grid._get_diagonal_list()
            self._node_nbrs = np.concatenate((self._node_nbrs, diag_nbrs), 1)
            self._link_lengths = np.empty(8, dtype=float)
            self._link_lengths[0] = dx
            self._link_lengths[2] = dx
            self._link_lengths[1] = dy
            self._link_lengths[3] = dy
            self._link_lengths[4:].fill(np.sqrt(dx*dx + dy*dy))
        elif ((type(self._grid) is landlab.grid.raster.RasterModelGrid) and
                (self._routing is 'D4')):
            self._link_lengths = np.empty(4, dtype=float)
            self._link_lengths[0] = dx
            self._link_lengths[2] = dx
            self._link_lengths[1] = dy
            self._link_lengths[3] = dy
        else:
            self._link_lengths = self._grid._length_of_link_with_diagonals
        self._lake_outlets = []  # a list of each unique lake outlet
        # ^note this is nlakes-long

        self.is_pit = self._grid.add_ones('node', 'is_pit', dtype=bool,
                                          noclobber=False)
        self.flood_status = self._grid.add_zeros('node', 'flood_status_code',
                                                 dtype=int, noclobber=False)
        self._lake_map = np.empty(self._grid.number_of_nodes, dtype=int)
        self._lake_map.fill(LOCAL_BAD_INDEX_VALUE)

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

        if type(self._grid) is landlab.grid.raster.RasterModelGrid:
            if not self._grid._diagonal_links_created:
                self._grid._setup_diagonal_links()

        h_diag = self._grid._diag_activelink_tonode
        t_diag = self._grid._diag_activelink_fromnode

        # These two lines assign the False flag to any node that is higher
        # than its partner on the other end of its link
        self.is_pit[h_orth[np.where(
            self._elev[h_orth] > self._elev[t_orth])[0]]] = False
        self.is_pit[t_orth[np.where(
            self._elev[t_orth] > self._elev[h_orth])[0]]] = False

        # If we have a raster grid, handle the diagonal active links too
        # (At the moment, their data structure is a bit different)
        # TODO: update the diagonal link data structures
        # DEJH doesn't understand why this can't be vectorized as above...
        if self._D8:
            for i in range(len(self._grid._diag_active_links)):
                h = self._grid._diag_activelink_tonode[i]
                t = self._grid._diag_activelink_fromnode[i]
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
                    elif self.flood_status[nbr] == _PIT or \
                            self.flood_status[nbr] == _FLOODED:
                        nodes_this_depression.append(nbr)
                        self.flood_status[nbr] = _CURRENT_LAKE
        return lowest_node

    def node_can_drain(self, the_node, nodes_this_depression):
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
        for nbr in self._node_nbrs[the_node]:
            if nbr != LOCAL_BAD_INDEX_VALUE:
                if self._elev[nbr] < self._elev[the_node] and \
                        self.flood_status[nbr] != _CURRENT_LAKE and \
                        self.flood_status[nbr] != _FLOODED:
                    # caveat about outlet elevation...
                    return True
        return False

    def is_valid_outlet(self, the_node, nodes_this_depression):
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

        if self.node_can_drain(the_node, nodes_this_depression):
            return True

    def _record_depression_depth_and_outlet(self, nodes_this_depression,
                                            outlet_id, pit_node):
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
            pit_node_where = np.searchsorted(self.pit_node_ids,
                                             pit_node)
            self._unique_pits[pit_node_where] = True
        elif np.any(fresh_nodes):  # lake is bigger than one or more existing
            self.flood_status[n] = _FLOODED
            depth_this_lake = self._elev[outlet_id] - self._elev[n]
            self.depression_depth[n] = depth_this_lake
            self.depression_outlet_map[n] = outlet_id
            # ^these two will just get stamped over as needed
            subsumed_lakes = np.unique(self._lake_map[n])  # IDed by pit_node
            # the final entry is LOCAL_BAD_INDEX_VALUE
            subs_lakes_where = np.searchsorted(self.pit_node_ids,
                                               subsumed_lakes[1:])
            pit_node_where = np.searchsorted(self.pit_node_ids,
                                             pit_node)
            self._unique_pits[subs_lakes_where] = False
            self._unique_pits[pit_node_where] = True
            self._pits_flooded -= (subsumed_lakes.size - 2)
            # -1 for the LOCAL_BAD_INDEX_VALUE that must be present; another -1
            # because a single lake is just replaced by a new lake.
            self._lake_map[n] = pit_node
        else:  # lake is subsumed within an existing lake
            assert np.all(np.equal(self.flood_status[n], _CURRENT_LAKE))
            self.flood_status[n] = _FLOODED
            pass

    def find_depression_from_pit(self, pit_node):
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
            lowest_node_on_perimeter = \
                self.find_lowest_node_on_lake_perimeter(nodes_this_depression)
            # note this can return the supplied node, if - somehow - the
            # surrounding nodes are all LOCAL_BAD_INDEX_VALUE
            found_outlet = self.is_valid_outlet(lowest_node_on_perimeter,
                                                nodes_this_depression)
            if not found_outlet:
                # Add lowest_node to the lake list
                nodes_this_depression.append(lowest_node_on_perimeter)
                # Flag it as being part of the current lake/depression
                self.flood_status[lowest_node_on_perimeter] = _CURRENT_LAKE
            # Safety check, in case a bug (ha!) puts us in an infinite loop
            assert (count < max_count), 'too many iterations in lake filler!'
            count += 1

        self.depression_outlets.append(lowest_node_on_perimeter)
        # Now that we've mapped this depression, record it in the arrays
        # depression_depth, depression_outlet, and flood_status
        self._record_depression_depth_and_outlet(nodes_this_depression,
                                                 lowest_node_on_perimeter,
                                                 pit_node)

        # TODO: ideally we need a way to keep track of the number, area extent,
        # and average depth of depressions. Tricky thing is that one might be
        # devoured by another, so would need to be removed from the list.

    def _identify_depressions_and_outlets(self):
        """Find depression and lakes on a topographic surface.

        Find and map the depressions/lakes in a topographic surface,
        given a previously identified list of pits (if any) in the surface.
        """
        self._pits_flooded = 0
        self._unique_pits = np.zeros_like(self.pit_node_ids, dtype=bool)
        for pit_node in self.pit_node_ids:
            self.find_depression_from_pit(pit_node)
            self._pits_flooded += 1
        assert len(self.depression_outlets) == self._unique_pits.size

        self.unique_lake_outlets = np.array(self.depression_outlets
                                            )[self._unique_pits]

    def map_depressions(self, pits='flow__sink_flag', reroute_flow=True):
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

        >>> rg = RasterModelGrid(5, 5)
        >>> z = rg.add_zeros('node', 'topographic__elevation')
        >>> z[:] = np.array([100., 100.,  95., 100., 100.,
        ...                  100., 101.,  92.,   1., 100.,
        ...                  100., 101.,   2., 101., 100.,
        ...                  100.,   3., 101., 101., 100.,
        ...                   90.,  95., 100., 100., 100.])
        >>> df = DepressionFinderAndRouter(rg)
        >>> df.map_depressions(pits=None, reroute_flow=False)
        >>> df.display_depression_map()
        . . . . .
        . . . ~ .
        . . ~ . .
        . ~ . . .
        o . . . .
        """
        self._lake_map.fill(LOCAL_BAD_INDEX_VALUE)
        self.depression_outlets = []  # reset these
        # Locate nodes with pits
        if type(pits) == str:
            try:
                pits = self._grid.at_node[pits]
                supplied_pits = np.where(pits)[0]
                self.pit_node_ids = as_id_array(
                    np.setdiff1d(supplied_pits, self._grid.boundary_nodes))
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
                np.setdiff1d(supplied_pits, self._grid.boundary_nodes))

            self.number_of_pits = self.pit_node_ids.size
            self.is_pit.fill(False)
            self.is_pit[self.pit_node_ids] = True
        # Set up "lake code" array
        self.flood_status.fill(_UNFLOODED)
        self.flood_status[self.pit_node_ids] = _PIT

        self._identify_depressions_and_outlets()

        if reroute_flow and ('flow__receiver_node' in
                             self._grid.at_node.keys()):
            self.receivers = self._grid.at_node['flow__receiver_node']
            self.sinks = self._grid.at_node['flow__sink_flag']
            self.grads = self._grid.at_node['topographic__steepest_slope']
            self._route_flow()
            self._reaccumulate_flow()

    def _route_flow(self):
        """Route flow across lake flats.

        Route flow across lake flats, which have already been identified.
        """
        for outlet_node, lake_code in zip(self.lake_outlets, self.lake_codes):
            nodes_in_lake = np.where(self.lake_map ==
                                     lake_code)[0]
            if len(nodes_in_lake) > 0:
                nodes_routed = np.array([outlet_node])
                # ^using set on assumption of cythonizing later
                nodes_on_front = np.array([outlet_node])
                self._handle_outlet_node(outlet_node, nodes_in_lake)
                while (len(nodes_in_lake) + 1) != len(nodes_routed):
                    if self._D8:
                        all_nbrs = np.hstack(
                            (self._grid.active_neighbors_at_node(
                             nodes_on_front),
                             self._grid._get_diagonal_list(
                             nodes_on_front)))
                    else:
                        all_nbrs = self._grid.active_neighbors_at_node(
                            nodes_on_front)
                    outlake = np.logical_not(np.in1d(all_nbrs.flat,
                                                     nodes_in_lake))
                    all_nbrs[outlake.reshape(all_nbrs.shape)] = -1
                    backflow = np.in1d(all_nbrs, nodes_routed)
                    all_nbrs[backflow.reshape(all_nbrs.shape)] = -1
                    (drains_from, unique_indxs) = np.unique(all_nbrs,
                                                            return_index=True)
                    # ^gets flattened, but, usefully, unique_indxs are
                    # *in order*
                    # remember, 1st one is always -1
                    # I bet this is sloooooooooow
                    good_nbrs = drains_from != -1
                    drains_to = nodes_on_front[unique_indxs[good_nbrs] //
                                               self.num_nbrs]
                    # to run the accumulator successfully, we need receivers,
                    # and sinks only. So the priority is sorting out the
                    # receiver field, and sealing the filled sinks (once while
                    # loop is done)
                    self.receivers[drains_from[good_nbrs]] = drains_to
                    # now put the relevant nodes in the relevant places:
                    nodes_on_front = drains_from[good_nbrs]
                    nodes_routed = np.union1d(nodes_routed, nodes_on_front)
                    self.grads[drains_from[good_nbrs]] = 0.
                    # ^downstream grad is 0.
        self.sinks[self.pit_node_ids] = False

    def _reaccumulate_flow(self):
        """Update drainage area, discharge, and upstream order.

        Invoke the accumulator a second time to update drainage area,
        discharge, and upstream order.
        """
        # Calculate drainage area, discharge, and downstr->upstr order
        Q_in = self._grid.at_node['water__unit_flux_in']
        areas = self._grid.cell_area_at_node.copy()
        areas[self._grid.closed_boundary_nodes] = 0.
        self.a, q, s = flow_accum_bw.flow_accumulation(self.receivers,
                                                       np.where(self.sinks)[0],
                                                       node_cell_area=areas,
                                                       runoff_rate=Q_in)
        # finish the property updating:
        self._grid.at_node['drainage_area'][:] = self.a
        self._grid.at_node['water__discharge'][:] = q
        self._grid.at_node['flow__upstream_node_order'][:] = s
        # ## TODO: No obvious easy way to recover the receiver_link.
        # ## Think more on this.
        # ## Right now, we're just not updating it.

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
                    (self._grid.active_neighbors_at_node(
                     outlet_node, bad_index=-1),
                     self._grid._get_diagonal_list(
                     outlet_node, bad_index=-1)))
            else:
                outlet_neighbors = self._grid.active_neighbors_at_node(
                    outlet_node, bad_index=-1).copy()
            inlake = np.in1d(outlet_neighbors.flat, nodes_in_lake)
            assert inlake.size > 0
            outlet_neighbors[inlake] = -1
            unique_outs, unique_indxs = np.unique(outlet_neighbors,
                                                  return_index=True)
            out_draining = unique_outs[1:]
            if isinstance(self._grid, landlab.grid.raster.RasterModelGrid):
                link_l = self._link_lengths
            else:  # Voronoi
                link_l = self._link_lengths[
                    self._grid.links_at_node[outlet_node, :]]
            eff_slopes = ((self._elev[outlet_node] -
                           self._elev[out_draining]) /
                          link_l[unique_indxs[1:]])
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
                    print('o', end=' ')
                elif self.flood_status[n] == _UNFLOODED:
                    print('.', end=' ')
                else:
                    print('~', end=' ')
                n += 1
            print()

    @property
    def lake_outlets(self):
        """
        Returns the *unique* outlets for each lake, in same order as the
        return from lake_codes.
        """
        return np.array(self.depression_outlets)[self._unique_pits]

    @property
    def lake_codes(self):
        """
        Returns the *unique* code assigned to each unique lake. These are
        the values used to map the lakes in the property "lake_map".
        """
        return self.pit_node_ids[self._unique_pits]

    @property
    def number_of_lakes(self):
        """
        Return the number of individual lakes.
        """
        return self._unique_pits.sum()

    @property
    def lake_map(self):
        """
        Return an array of ints, where each node within a lake is labelled
        with a unique (non-consecutive) code corresponding to each unique
        lake. The codes used can be obtained with *lake_codes*.
        Nodes not in a lake are labelled with LOCAL_BAD_INDEX_VALUE.
        """
        return self._lake_map

    @property
    def lake_at_node(self):
        """
        Return a boolean array, True if the node is flooded, False otherwise.
        """
        return self._lake_map != LOCAL_BAD_INDEX_VALUE

    @property
    def lake_areas(self):
        """
        A nlakes-long array of the area of each lake. The order is the same as
        that returned by *lake_codes*.
        """
        lake_areas = np.empty(self.number_of_lakes)
        lake_counter = 0
        for lake_code in self.lake_codes:
            each_cell_in_lake = self._grid.cell_area_at_node[self.lake_map ==
                                                             lake_code]
            lake_areas[lake_counter] = each_cell_in_lake.sum()
            lake_counter += 1
        return lake_areas

    @property
    def lake_volumes(self):
        """
        A nlakes-long array of the volume of each lake. The order is the same
        as that returned by *lake_codes*.
        """
        lake_vols = np.empty(self.number_of_lakes)
        lake_counter = 0
        col_vols = self._grid.cell_area_at_node * self.depression_depth
        for lake_code in self.lake_codes:
            each_cell_in_lake = col_vols[self.lake_map == lake_code]
            lake_vols[lake_counter] = each_cell_in_lake.sum()
            lake_counter += 1
        return lake_vols


def main():
    """temporary: test."""
    print('howdy')
    from landlab import RasterModelGrid
    from numpy.random import rand
    grid = RasterModelGrid(4, 5, 1.0)
    z = grid.add_zeros('node', 'topographic__elevation')
    z[:] = rand(len(z)) * 100
    print(z)
    dep_finder = DepressionFinderAndRouter(grid,  '/Users/gtucker/Dev/' +
                                           'Landlab/gt_tests/test_inputs_for' +
                                           '_depression_mapper.txt')
    dep_finder.map_depressions()

    n = 0
    for r in range(grid.number_of_node_rows - 1, -1, -1):
        for c in range(grid.number_of_node_columns):
            print(int(z[n]), '(', dep_finder.is_pit[n], ')',  end=' ')
            n += 1
        print()

    n = 0
    for r in range(grid.number_of_node_rows):
        for c in range(grid.number_of_node_columns):
            if dep_finder.depression_outlet_map[n] == LOCAL_BAD_INDEX_VALUE:
                print(dep_finder.depression_depth[n], '( X )', end=' ')
            else:
                print(dep_finder.depression_depth[n], '(',
                      dep_finder.depression_outlet_map[n], ')',  end=' ')
            n += 1
        print()

    dep_finder.display_depression_map()

    # Now, route flow through.
    # First, find flow dirs without lakes. Then, adjust.
    from landlab.components.flow_routing import grid_flow_directions
    (rcvr, ss) = grid_flow_directions(grid, z)
    print(z)
    print(rcvr)
    print(ss)


if __name__ == '__main__':
    import doctest
    doctest.testmod()
    main()

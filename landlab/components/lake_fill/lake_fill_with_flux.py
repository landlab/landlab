#!/usr/env/python

"""
lake_filler_with_flux.py: Component to fill depressions in a landscape while
honouring mass balance.

Similar to the DepressionFinderAndRouter, but will not fill a lake to the brim
if there is not enough incoming flux to do so. Designed to "play nice" with
the FlowAccumulator.
"""

from __future__ import print_function

import warnings

from landlab import FieldError, Component
from landlab import RasterModelGrid, VoronoiDelaunayGrid  # for type tests
from landlab.components import LakeMapperBarnes
from landlab.components.lake_fill import StablePriorityQueue
from landlab.utils.return_array import return_array_at_node
from landlab.core.messages import warning_message

from landlab import BAD_INDEX_VALUE
import six
import numpy as np

LOCAL_BAD_INDEX_VALUE = BAD_INDEX_VALUE
LARGE_ELEV = 9999999999.


def fill_depression_from_pit_discharges(mg, depression_outlet, depression_nodes,
                                        pit_nodes, Vw_at_pits, surface_z,
                                        water_z, water_vol_balance_terms,
                                        neighbors_at_nodes):
    """
    Take an outlet and its catchment nodes, then work outwards from the
    pit nodes using the volumes of water at each to work out which nodes
    are connected to which, and where is inundated.

    Assumes we already know that the pit doesn't overtop, i.e., is its own
    little endorheic basin.

    DO NOT INCLINE THE WATER SURFACES!!

    Idea here is that water_vol_balance_terms IS A FUNC OR ARRAY, such that
    it can be sensitive to depth if necessary.
    Provide water_vol_balance_terms as (c, K), where balance = K*depth + c.
    Either can be 0. Equation must be linear. Note both terms are likely to
    be negative... Much safer to ensure these always are, and add any
    positive terms as part of the accumulation algorithm inputs
    (AssertionErrors are likely for large gains...)
    """
    disch_map = mg.zeros('node', dtype=float)  # really a vol map!
    lake_map = mg.zeros('node', dtype=int)
    lake_map.fill(-1)  # this will store the unique lake IDs
    area_map = mg.cell_area_at_node  # we will update this as we go...
    # sort the pit nodes according to existing fill.
    # We need to start at the lowest.
    num_pits = len(pit_nodes)
    zsurf_pit_nodes = water_z[pit_nodes]  # used to be surface_z. Does it matter?
    pit_node_sort = np.argsort(zsurf_pit_nodes)
    zwater_nodes_in_order = water_z[pit_node_sort]
    pit_nodes_in_order = pit_nodes[pit_node_sort]
    Vw_at_pits_in_order = Vw_at_pits[pit_node_sort]
    vol_rem_at_pit = {}
    accum_area_at_pit = {}
    lake_water_level = {}
    accum_K_at_pit = {}
    lake_is_full = {}
    for pit, vol, A, level in zip(pit_nodes_in_order, Vw_at_pits_in_order,
                                  area_map[pit_node_sort],
                                  zwater_nodes_in_order):
        vol_rem_at_pit[pit] = vol
        lake_water_level[pit] = level
        lake_is_full[pit] = False
        accum_area_at_pit[pit] = 0.
        accum_K_at_pit[pit] = 0.


def raise_lakes_to_water_exhaustion(
        flow__receiver_nodes,
        link_lengths,
        links_at_node,
        lake_map,
        area_map,
        master_pit_q,
        lake_q_dict,
        init_water_surface,
        lake_water_level,
        accum_area_at_pit,
        vol_rem_at_pit,
        accum_K_at_pit,
        lake_is_full,
        lake_spill_node,
        water_vol_balance_terms,
        neighbors,
        closednodes,
        drainingnodes,
        water_out_dict):
    """
    For known water volumes at pits, raise a water surface to honour those
    volumes and any gains and losses created by the lakes themselves.

    Parameters
    ----------
    flow__receiver_nodes : array of ints
        The nodes already known to be the steepest descent path on the
        unflooded terrain. Typically the output from a many-to-one
        FlowDirector.
    link_lengths : array of floats
        The length of each link on the grid. Typically grid.length_of_link.
    links_at_node : (nnodes, max_nlinks_at_node)
        The links at each node on the grid, -1 where undefined. Typically
        grid.links_at_node or d8 equivalent.
    lake_map : array of ints
        Grid-node array of -1s with IDs where a node is known to be flooded
        from that pit, and ID-related coding for lake adjacent or lake outlet
        nodes.
    area_map : array of floats
        Grid-node array of the areas of each cell; typically the output of
        a grid cell_area_at_node property.
    master_pit_q: StablePriorityQueue
        This queue holds the current lakes in the grid, ordered by water
        surface elevation priority. The current lake will have just been
        popped off this list, exposing the next highest surface.
    lake_q_dict : dict of StablePriorityQueues
        A dict of queues holding nodes the node which begins the iteration
        (typically the pit node), along with the first set of nodes to
        explore. If the lake hasn't been filled before, this will contain
        simply the pit node, and its immediate neighbors.
    init_water_surface : array
        The original topo of the water surface at the start of the step,
        i.e., before we start moving the lake water levels this time around.
        This will follow the land surface outside any pre-existing lakes.
    lake_water_level : dict
        The current (transient) water level at each lake.
    accum_area_at_pit : dict
        The total surface area of each lake, updated through time.
    vol_rem_at_pit : dict
        The "excess" water volume available at each lake, that will be used to
        continue to raise the water level until exhausted.
    accum_K_at_pit : dict
        Stores the current total K term in the (optional) water balance
        equation across the sum of all nodes in each lake. Needed to perform
        water balance as the lake level rises during the fill.
    lake_is_full : dict
        Stores the status of each distinct pit in the depression.
    lake_spill_node : dict
        the current spill node of each lake. -1 if not defined.
    water_vol_balance_terms : (float c, float K)
        polyparams for a linear fit of the water vol balance func, such that
        dV = A * (K*z + c). Fit must be linear.
    neighbors : (nnodes, max_nneighbors) array
        The neighbors at each node on the grid. Typically obtained direct
        from a grid, as the adjacent_nodes_at_node or the equivalent
        d8 structure.
    closednodes : array of bool
        Nodes not to be explored by the algorithm.
    drainingnodes : array of bool
        Once a lake raises itself to the level of one of these nodes, it will
        cease to rise any more, ever (i.e., these are boundary nodes!)
    water_out_dict : dict
        Dict with the True drainingnodes as keys and the water volumes
        reaching those nodes as values. Only gets populated where lake
        overflow discharges actually reach the node; also note that this
        water volume is only what is overtopped, and it does not include
        any discharge accumulating at that node "conventionally".

    Examples
    --------
    >>> import numpy as np
    >>> from landlab.components.lake_fill import StablePriorityQueue
    >>> from landlab.components import FlowDirectorSteepest
    >>> from landlab import RasterModelGrid
    >>> init_water_surface = np.array([  0.,  0.,  0.,  0.,  0.,  0.,
    ...                                  0., -9., -6.,  0., -8.,  0.,
    ...                                  0., -8., -5., -6., -7.,  0.,
    ...                                  0., -7.,  0., -7.,  0.,  0.,
    ...                                  0.,  0.,  0., -8., -3., -1.,
    ...                                  0.,  0.,  0.,  0., -2.,  0.])
    >>> lake_map = -1 * np.ones(36, dtype=int)
    >>> area_map = 4. * np.ones(36, dtype=float)
    >>> master_pit_q = StablePriorityQueue()
    >>> lake_q_dict = {7: StablePriorityQueue(), 10: StablePriorityQueue(),
    ...                27: StablePriorityQueue()}
    >>> lake_water_level = {7: init_water_surface[7],
    ...                     10: init_water_surface[10],
    ...                     27: init_water_surface[27]}
    >>> accum_area_at_pit = {7: 0., 10: 0., 27: 0.}
    >>> vol_rem_at_pit = {7: 200., 10: 8., 27: 17.}
    >>> accum_K_at_pit = {7: 0., 10: 0., 27: 0.}
    >>> lake_is_full = {7: False, 10: False, 27: False}
    >>> lake_spill_node = {7: -1, 10: -1, 27: -1}
    >>> grid = RasterModelGrid((6, 6), 2.)
    >>> _ = grid.add_field(
    ...     'node', 'topographic__elevation', init_water_surface)
    >>> fd = FlowDirectorSteepest(grid)
    >>> fd.run_one_step()
    >>> flow__receiver_nodes = grid.at_node['flow__receiver_node']
    >>> link_lengths = grid.length_of_link
    >>> links_at_node = grid.links_at_node
    >>> neighbors = grid.adjacent_nodes_at_node
    >>> closednodes = np.zeros(36, dtype=bool)
    >>> closednodes[34] = True
    >>> drainingnodes = np.zeros(36, dtype=bool)
    >>> drainingnodes[29] = True  # create an outlet
    >>> water_vol_balance_terms = (1., -4.)
    >>> water_out_dict = {}

    >>> for cpit in (7, 10, 27):
    ...     lake_q_dict[cpit].add_task(cpit, priority=lake_water_level[cpit])

    >>> for cpit in (7, 10, 27):
    ...     for nghb in neighbors[cpit]:
    ...         lake_q_dict[cpit].add_task(nghb,
    ...                                    priority=init_water_surface[nghb])
    ...         lake_map[nghb] = cpit + 36

    >>> for cpit in (7, 10, 27):
    ...     elev = lake_water_level[cpit]
    ...     master_pit_q.add_task(cpit, priority=elev)

    >>> raise_lakes_to_water_exhaustion(
    ...     flow__receiver_nodes, link_lengths, links_at_node, lake_map,
    ...     area_map, master_pit_q, lake_q_dict, init_water_surface,
    ...     lake_water_level, accum_area_at_pit, vol_rem_at_pit,
    ...     accum_K_at_pit, lake_is_full, lake_spill_node,
    ...     water_vol_balance_terms, neighbors, closednodes, drainingnodes,
    ...     water_out_dict)

    >>> np.all(np.equal(lake_map, np.array([-1, 43, 43, -1, 46, -1,
    ...                                     43,  7,  7, 46, 10, 46,
    ...                                     43,  7, 79, 63, 10, 46,
    ...                                     43,  7, 63, 27, 63, -1,
    ...                                     -1, 43, 63, 27, 63, -1,
    ...                                     -1, -1, -1, 63, -1, -1])))
    True
    >>> len(master_pit_q.tasks_currently_in_queue())  # we're all done
    0
    >>> lake_q_dict[7].tasks_currently_in_queue()  # kept just in case
    array([ 6, 12,  1, 18,  2, 20, 25,  9])
    >>> np.allclose(np.sort(lake_water_level.values()), (-6.875, -6.7, -5))
    True
    >>> accum_area_at_pit == {7: 16.0, 10: 8.0, 27: 8.0}
    True
    >>> vol_rem_at_pit == {7: 0.0, 10: 0.0, 27: 0.0}
    True
    >>> accum_K_at_pit == {7: -64.0, 10: -16.0, 27: -16.0}
    True
    >>> lake_is_full == {7: False, 10: True, 27: True}
    True
    >>> lake_spill_node == {7: 14, 10: -1, 27: -1}
    True
    >>> water_out_dict == {}  # no water is escaping
    True

    Note here that we've engineered a case where even though lake 7 is
    overtopping, the additional discharge is still not enough to raise the
    downstream lake to its level, so we retain all three distinct lakes.

    Now juice the system up again:

    >>> vol_rem_at_pit[7] += 2000.  # recharge
    >>> master_pit_q.add_task(7)  # re-permit action on the lake
    >>> raise_lakes_to_water_exhaustion(
    ...     flow__receiver_nodes, link_lengths, links_at_node, lake_map,
    ...     area_map, master_pit_q, lake_q_dict, init_water_surface,
    ...     lake_water_level, accum_area_at_pit, vol_rem_at_pit,
    ...     accum_K_at_pit, lake_is_full, lake_spill_node,
    ...     water_vol_balance_terms, neighbors, closednodes, drainingnodes,
    ...     water_out_dict)

    >>> np.all(np.equal(lake_map, np.array([-1, 43, 43, -1, 43, -1,
    ...                                     43,  7,  7, 43,  7, 43,
    ...                                     43,  7,  7,  7,  7, 43,
    ...                                     43,  7, 43,  7, 43, -1,
    ...                                     -1, 43, 43,  7,  7, 79,
    ...                                     -1, -1, -1, 43, -1, -1])))
    True
    >>> len(master_pit_q.tasks_currently_in_queue())  # we're all done
    0
    >>> lake_q_dict[7].tasks_currently_in_queue()  # kept just in case
    array([11,  6, 33,  4,  9, 26,  1,  2, 22, 20, 12, 17, 25, 18])
    >>> lake_water_level == {7: -1.0}
    True
    >>> accum_area_at_pit == {7: 44.0}
    True
    >>> vol_rem_at_pit == {7: 1009.0}
    True
    >>> accum_K_at_pit == {7: -176.0}
    True
    >>> lake_is_full == {7: False}
    True
    >>> lake_spill_node == {7: 29}
    True
    >>> water_out_dict == {29: 1009.0}  # note the same as vol_rem_at_pit
    True

    In this second example, we demonstrate that this function can correctly
    cascade water down a chain of lakes:

    >>> from landlab import CLOSED_BOUNDARY
    >>> init_water_surface = np.array(
    ...     [  0.,  0.,  0.,  0.,  0.,  0.,  0.,  0.,  0.,
    ...        0.,  9., -2.,  8.,  4.,  5., -1.,  6.,  0.,
    ...        0.,  0.,  0.,  0.,  0.,  0.,  0.,  0.,  0.])
    >>> lake_map = -1 * np.ones(27, dtype=int)
    >>> area_map = 4. * np.ones(27, dtype=float)
    >>> master_pit_q = StablePriorityQueue()
    >>> mypits = (11, 13, 15)
    >>> lake_q_dict = {}
    >>> lake_water_level = {}
    >>> accum_area_at_pit = {}
    >>> accum_K_at_pit = {}
    >>> lake_is_full = {}
    >>> lake_spill_node = {}
    >>> for cpit in mypits:
    ...     lake_water_level[cpit] = init_water_surface[cpit]
    ...     lake_q_dict[cpit] = StablePriorityQueue()
    ...     accum_area_at_pit[cpit] = 0.
    ...     accum_K_at_pit[cpit] = 0.
    ...     lake_is_full[cpit] = False
    ...     lake_spill_node[cpit] = -1
    >>> vol_rem_at_pit = {11: 100., 13: 100., 15: 100.}
    >>> grid = RasterModelGrid((3, 9), 2.)
    >>> _ = grid.add_field(
    ...     'node', 'topographic__elevation', init_water_surface)
    >>> grid.status_at_node[:9] = CLOSED_BOUNDARY
    >>> grid.status_at_node[18:] = CLOSED_BOUNDARY
    >>> fd = FlowDirectorSteepest(grid)
    >>> fd.run_one_step()
    >>> flow__receiver_nodes = grid.at_node['flow__receiver_node']
    >>> link_lengths = grid.length_of_link
    >>> links_at_node = grid.links_at_node
    >>> neighbors = grid.adjacent_nodes_at_node
    >>> closednodes = np.zeros(27, dtype=bool)
    >>> closednodes[:9] = True
    >>> closednodes[18:] = True
    >>> drainingnodes = np.zeros(27, dtype=bool)
    >>> drainingnodes[[9, 17]] = True  # create competing outlets
    >>> water_vol_balance_terms = (0., 0.)
    >>> water_out_dict = {}

    >>> for cpit in mypits:
    ...     lake_q_dict[cpit].add_task(cpit, priority=lake_water_level[cpit])

    >>> for cpit in mypits:
    ...     for nghb in neighbors[cpit]:
    ...         if not closednodes[nghb]:
    ...             lake_q_dict[cpit].add_task(
    ...                 nghb, priority=init_water_surface[nghb])
    ...             lake_map[nghb] = cpit + 27

    >>> for cpit in mypits:
    ...     elev = lake_water_level[cpit]
    ...     master_pit_q.add_task(cpit, priority=elev)

    >>> raise_lakes_to_water_exhaustion(
    ...     flow__receiver_nodes, link_lengths, links_at_node, lake_map,
    ...     area_map, master_pit_q, lake_q_dict, init_water_surface,
    ...     lake_water_level, accum_area_at_pit, vol_rem_at_pit,
    ...     accum_K_at_pit, lake_is_full, lake_spill_node,
    ...     water_vol_balance_terms, neighbors, closednodes, drainingnodes,
    ...     water_out_dict)

    >>> np.all(np.equal(lake_map, np.array(
    ...     [-1, -1, -1, -1, -1, -1, -1, -1, -1,
    ...      -1, 38, 11, 65, 15, 15, 15, 69, -1,
    ...      -1, -1, -1, -1, -1, -1, -1, -1, -1])))
    True
    >>> lake_water_level == {11: 8.0, 15: 6.0}
    True
    >>> accum_area_at_pit == {11: 4.0, 15: 12.0}
    True
    >>> vol_rem_at_pit == {11: 0.0, 15: 0.0}
    True
    >>> lake_is_full == {11: False, 15: False}
    True
    >>> lake_spill_node == {11: 12, 15: 16}
    True
    >>> water_out_dict == {17: 220.0}
    True
    """
    flow_goes_to = {}  # this is for cascading flow
    while True:
        # print('outer loop')
        try:
            cpit = master_pit_q.pop_task()  # check the lowest node
        except KeyError:  # all done
            break
        # Call below will break when water at that pit is exhausted, or
        # when it reaches a sill. In the latter case, the pit will have
        # been put back in the master queue.
        # These are the important cases. If one of these gets popped off
        # the queue here, it means that the *lowest* thing was a sill.
        # That means rather than raise it, we need to flood the water down
        # and out of it; i.e., we should NEVER be raising a lake at a sill.
        # So before we call _raise_lake_to_limit, first check the
        # lake_spill_node:
        if cpit not in lake_spill_node.keys():
            # possible to pop a lake that doesn't exist any more, so
            # print('No more lake!')
            continue
        if lake_spill_node[cpit] != -1:
            # print('spilling')
            # flood water down
            next_lake, escape_ID = _route_outlet_to_next(
                lake_spill_node[cpit], cpit, flow__receiver_nodes,
                link_lengths, links_at_node, init_water_surface, lake_map,
                neighbors, closednodes, drainingnodes)
            flow_goes_to[cpit] = (next_lake, escape_ID)
            # print('nextlake:', next_lake)
            while next_lake != -1:
                next_lake = lake_map[next_lake]
                # ^^in case lakes have merged since flow_goes_to was updated
                # this will now be -1 if we're draining to a drainingnode,
                # as well as if escape_ID was returned above
                try:
                    next_lake, escape_ID = flow_goes_to[next_lake]
                except KeyError:
                    break
                    # we don't know where next_lake drains to, or it just
                    # doesn't drain yet
            if next_lake != -1:
                vol_rem_at_pit[next_lake] = vol_rem_at_pit[cpit]
                vol_rem_at_pit[cpit] = 0.
                lake_is_full[next_lake] = False
                # now operate immediately on this node
                cpit = next_lake
            else:
                # we are draining out. If we get here, there's no way this
                # lake is of any importance any more to the lakes we are
                # still working on, so
                try:
                    water_out_dict[escape_ID] += vol_rem_at_pit[cpit]
                    # this necessary in case there is already water spilling,
                    # and more is added later from above (a cascade)
                except KeyError:
                    water_out_dict[escape_ID] = vol_rem_at_pit[cpit]
                vol_rem_at_pit[cpit] = 0.
                continue
        else:
            pass

        _raise_lake_to_limit(cpit, lake_map, area_map,
                             master_pit_q,
                             lake_q_dict, init_water_surface,
                             lake_water_level,
                             accum_area_at_pit,
                             vol_rem_at_pit,
                             accum_K_at_pit,
                             lake_is_full,
                             lake_spill_node,
                             water_vol_balance_terms,
                             neighbors, closednodes, drainingnodes,
                             water_out_dict)


def _raise_lake_to_limit(current_pit, lake_map, area_map,
                         master_pit_q,
                         lake_q_dict, init_water_surface,
                         lake_water_level,
                         accum_area_at_pit,
                         vol_rem_at_pit,
                         accum_K_at_pit,
                         lake_is_full,
                         lake_spill_node,
                         water_vol_balance_terms,
                         neighbors, closednodes, drainingnodes,
                         water_out_dict):
    """
    Lift a lake level from a starting elevation to a new break point.
    Break points are specified by (a) running out of discharge, (b)
    making contact with another lake,
    or (c) reaching a spill point for the current depression. (b) is a bit
    of a special case - we test for this by contact, not by level. We are
    permitted to keep raising the level past the lake level if we don't
    actually touch another lake, don't meet a spill, & still have enough
    water to do it.

    Water balance losses are handled in a fairly unsophisticated manner,
    largely for speed. Large losses may result in violations of water balance.

    Parameters
    ----------
    Current_pit : int
        Node ID of the pit that uniquely identifies the current lake that
        we are filling.
    lake_map : array of ints
        Grid-node array of -1s with IDs where a node is known to be flooded
        from that pit, and ID-related coding for lake adjacent or lake outlet
        nodes.
    area_map : array of floats
        Grid-node array of the areas of each cell; typically the output of
        a grid cell_area_at_node property.
    master_pit_q: StablePriorityQueue
        This queue holds the current lakes in the grid, ordered by water
        surface elevation priority. The current lake will have just been
        popped off this list, exposing the next highest surface.
    lake_q_dict : dict of StablePriorityQueues
        A dict of queues holding nodes the node which begins the iteration
        (typically the pit node), along with the first set of nodes to
        explore. We hold this information through time since we return to
        each lake repeatedly to update the elevs. If the lake hasn't been
        filled before, this will contain simply the pit node, and its
        immediate neighbors.
    init_water_surface : array
        The original topo of the water surface at the start of the step,
        i.e., before we start moving the lake water levels this time around.
        This will follow the land surface outside any pre-existing lakes.
    lake_water_level : dict
        The current (transient) water level at each lake.
    accum_area_at_pit : dict
        The total surface area of each lake, updated through time.
    vol_rem_at_pit : dict
        The "excess" water volume available at each lake, that will be used to
        continue to raise the water level until exhausted.
    accum_K_at_pit : dict
        Stores the current total K term in the (optional) water balance
        equation across the sum of all nodes in each lake. Needed to perform
        water balance as the lake level rises during the fill.
    lake_is_full : dict
        Stores the status of each distinct pit in the depression.
    lake_spill_node : dict
        the current spill node of each lake. -1 if not defined.
    water_vol_balance_terms : (float c, float K)
        polyparams for a linear fit of the water vol balance func, such that
        dV = A * (K*z + c). Fit must be linear.
    neighbors : (nnodes, max_nneighbors) array
        The neighbors at each node on the grid. Typically obtained direct
        from a grid, as the adjacent_nodes_at_node or the equivalent
        d8 structure.
    closednodes : array of bool
        Nodes not to be explored by the algorithm.
    drainingnodes : array of bool
        Once a lake raises itself to the level of one of these nodes, it will
        cease to rise any more, ever (i.e., these are boundary nodes!)
    water_out_dict : dict
        Dict with the True drainingnodes as keys and the water volumes
        reaching those nodes as values. Only gets populated where lake
        overflow discharges actually reach the node; also note that this
        water volume is only what is overtopped, and it does not include
        any discharge accumulating at that node "conventionally".

    Examples
    --------
    >>> import numpy as np
    >>> from landlab.components.lake_fill import StablePriorityQueue
    >>> from landlab import RasterModelGrid
    >>> init_water_surface = np.array([  0.,  0.,  0.,  0.,  0.,  0.,
    ...                                  0., -9., -6.,  0., -8.,  0.,
    ...                                  0., -8., -5., -6., -7.,  0.,
    ...                                  0., -7.,  0., -7.,  0.,  0.,
    ...                                  0.,  0.,  0., -8., -3., -1.,
    ...                                  0.,  0.,  0.,  0., -2.,  0.])
    >>> lake_map = -1 * np.ones(36, dtype=int)
    >>> area_map = 4. * np.ones(36, dtype=float)
    >>> master_pit_q = StablePriorityQueue()
    >>> lake_q_dict = {7: StablePriorityQueue(), 10: StablePriorityQueue(),
    ...                27: StablePriorityQueue()}
    >>> lake_water_level = {7: init_water_surface[7],
    ...                     10: init_water_surface[10],
    ...                     27: init_water_surface[27]}
    >>> accum_area_at_pit = {7: 0., 10: 0., 27: 0.}
    >>> vol_rem_at_pit = {7: 200., 10: 8., 27: 17.}
    >>> accum_K_at_pit = {7: 0., 10: 0., 27: 0.}
    >>> lake_is_full = {7: False, 10: False, 27: False}
    >>> lake_spill_node = {7: -1, 10: -1, 27: -1}
    >>> grid = RasterModelGrid((6, 6), 2.)
    >>> neighbors = grid.adjacent_nodes_at_node
    >>> closednodes = np.zeros(36, dtype=bool)
    >>> closednodes[34] = True
    >>> drainingnodes = np.zeros(36, dtype=bool)
    >>> drainingnodes[29] = True  # create an outlet
    >>> water_vol_balance_terms = (0., 0.)
    >>> water_out_dict = {}

    Now get the lake initial conditions sorted out:

    >>> for cpit in (7, 10, 27):
    ...     lake_q_dict[cpit].add_task(cpit, priority=lake_water_level[cpit])

    >>> for cpit in (7, ):
    ...     for nghb in neighbors[cpit]:
    ...         lake_q_dict[cpit].add_task(nghb,
    ...                                    priority=init_water_surface[nghb])
    ...         lake_map[nghb] = cpit + 36

    >>> _raise_lake_to_limit(7, lake_map, area_map, master_pit_q,
    ...                      lake_q_dict, init_water_surface, lake_water_level,
    ...                      accum_area_at_pit, vol_rem_at_pit,
    ...                      accum_K_at_pit, lake_is_full, lake_spill_node,
    ...                      water_vol_balance_terms, neighbors, closednodes,
    ...                      drainingnodes, water_out_dict)

    >>> np.all(np.equal(lake_map, np.array([-1, 43, 43, -1, -1, -1,
    ...                                     43,  7,  7, 43, -1, -1,
    ...                                     43,  7, 79, -1, -1, -1,
    ...                                     43,  7, 43, -1, -1, -1,
    ...                                     -1, 43, -1, -1, -1, -1,
    ...                                     -1, -1, -1, -1, -1, -1])))
    True
    >>> master_pit_q.tasks_currently_in_queue()  # sill, so 7 is in it after
    array([7])
    >>> lake_q_dict[7].tasks_currently_in_queue()
    array([ 6, 12,  1, 18,  2, 20, 25,  9])
    >>> lake_water_level == {7: -5.0, 10: -8.0, 27: -8.0}
    True
    >>> accum_area_at_pit == {7: 16.0, 10: 0.0, 27: 0.0}
    True
    >>> vol_rem_at_pit == {7: 160.0, 10: 8.0, 27: 17.0}
    True
    >>> accum_K_at_pit == {7: 0.0, 10: 0.0, 27: 0.0}
    True
    >>> lake_is_full == {7: False, 10: False, 27: False}
    True
    >>> lake_spill_node == {7: 14, 10: -1, 27: -1}
    True

    >>> for cpit in (10, ):
    ...     for nghb in neighbors[cpit]:
    ...         lake_q_dict[cpit].add_task(nghb,
    ...                                    priority=init_water_surface[nghb])
    ...         lake_map[nghb] = cpit + 36

    >>> _raise_lake_to_limit(10, lake_map, area_map, master_pit_q,
    ...                      lake_q_dict, init_water_surface, lake_water_level,
    ...                      accum_area_at_pit, vol_rem_at_pit,
    ...                      accum_K_at_pit, lake_is_full, lake_spill_node,
    ...                      water_vol_balance_terms, neighbors, closednodes,
    ...                      drainingnodes, water_out_dict)

    >>> np.all(np.equal(lake_map, np.array([-1, 43, 43, -1, 46, -1,
    ...                                     43,  7,  7, 46, 10, 46,
    ...                                     43,  7, 79, 46, 10, 46,
    ...                                     43,  7, 43, -1, 46, -1,
    ...                                     -1, 43, -1, -1, -1, -1,
    ...                                     -1, -1, -1, -1, -1, -1])))
    True
    >>> master_pit_q.tasks_currently_in_queue() # 10 filled, so left off
    array([7])
    >>> lake_q_dict[10].tasks_currently_in_queue()  # half-filled node back in
    array([16,  9, 15, 17, 22,  4, 11])
    >>> lake_water_level == {7: -5.0, 10: -6.5, 27: -8.0}
    True
    >>> accum_area_at_pit == {7: 16.0, 10: 8.0, 27: 0.0}
    True
    >>> vol_rem_at_pit == {7: 160.0, 10: 0.0, 27: 17.0}
    True
    >>> accum_K_at_pit == {7: 0.0, 10: 0.0, 27: 0.0}
    True
    >>> lake_is_full == {7: False, 10: True, 27: False}
    True
    >>> lake_spill_node == {7: 14, 10: -1, 27: -1}
    True

    >>> vol_rem_at_pit[10] = 9.  # recharge to continue fill
    >>> _raise_lake_to_limit(10, lake_map, area_map, master_pit_q,
    ...                      lake_q_dict, init_water_surface, lake_water_level,
    ...                      accum_area_at_pit, vol_rem_at_pit,
    ...                      accum_K_at_pit, lake_is_full, lake_spill_node,
    ...                      water_vol_balance_terms, neighbors, closednodes,
    ...                      drainingnodes, water_out_dict)

    >>> np.all(np.equal(lake_map, np.array([-1, 43, 43, -1, 46, -1,
    ...                                     43,  7,  7, 46, 10, 46,
    ...                                     43,  7, 79, 82, 10, 46,
    ...                                     43,  7, 43, -1, 46, -1,
    ...                                     -1, 43, -1, -1, -1, -1,
    ...                                     -1, -1, -1, -1, -1, -1])))
    True
    >>> master_pit_q.tasks_currently_in_queue()
    array([10, 7])
    >>> lake_q_dict[10].tasks_currently_in_queue()
    array([14,  9, 11, 17, 22,  4])
    >>> lake_water_level == {7: -5.0, 10: -6.0, 27: -8.0}
    True
    >>> accum_area_at_pit == {7: 16.0, 10: 8.0, 27: 0.0}
    True
    >>> vol_rem_at_pit == {7: 160.0, 10: 5.0, 27: 17.0}
    True
    >>> lake_is_full == {7: False, 10: False, 27: False}
    True
    >>> lake_spill_node == {7: 14, 10: 15, 27: -1}
    True

    Now, if we call this again, nothing happens at all, as we're at a sill
    already:

    >>> _raise_lake_to_limit(10, lake_map, area_map, master_pit_q,
    ...                      lake_q_dict, init_water_surface, lake_water_level,
    ...                      accum_area_at_pit, vol_rem_at_pit,
    ...                      accum_K_at_pit, lake_is_full, lake_spill_node,
    ...                      water_vol_balance_terms, neighbors, closednodes,
    ...                      drainingnodes, water_out_dict)

    >>> np.all(np.equal(lake_map, np.array([-1, 43, 43, -1, 46, -1,
    ...                                     43,  7,  7, 46, 10, 46,
    ...                                     43,  7, 79, 82, 10, 46,
    ...                                     43,  7, 43, -1, 46, -1,
    ...                                     -1, 43, -1, -1, -1, -1,
    ...                                     -1, -1, -1, -1, -1, -1])))
    True
    >>> lake_water_level == {7: -5.0, 10: -6.0, 27: -8.0}
    True
    >>> lake_is_full == {7: False, 10: False, 27: False}
    True

    >>> for cpit in (27, ):
    ...     for nghb in neighbors[cpit]:
    ...         lake_q_dict[cpit].add_task(nghb,
    ...                                    priority=init_water_surface[nghb])
    ...         lake_map[nghb] = max((cpit + 36, lake_map[nghb]))
    >>> _raise_lake_to_limit(27, lake_map, area_map, master_pit_q,
    ...                      lake_q_dict, init_water_surface, lake_water_level,
    ...                      accum_area_at_pit, vol_rem_at_pit,
    ...                      accum_K_at_pit, lake_is_full, lake_spill_node,
    ...                      water_vol_balance_terms, neighbors, closednodes,
    ...                      drainingnodes, water_out_dict)

    The above call seamlessly merges the lakes & continues the subsequent rise
    until the (combined) water available craps out:

    >>> np.all(np.equal(lake_map, np.array([-1, 43, 43, -1, 46, -1,
    ...                                     43,  7,  7, 46, 10, 46,
    ...                                     43,  7, 79, 10, 10, 46,
    ...                                     43,  7, 46, 10, 46, -1,
    ...                                     -1, 43, 46, 10, 46, -1,
    ...                                     -1, -1, -1, 46, -1, -1])))
    True
    >>> lake_q_dict[10].tasks_currently_in_queue()
    array([15, 14, 11, 17, 28, 26,  4, 20, 22,  9, 33])
    >>> lake_water_level == {7: -5.0, 10: -5.5}
    True
    >>> accum_area_at_pit == {7: 16.0, 10: 20.0}
    True
    >>> vol_rem_at_pit == {7: 160.0, 10: 0.0}
    True
    >>> lake_is_full == {7: False, 10: True}
    True
    >>> lake_spill_node == {7: 14, 10: -1}
    True

    >>> vol_rem_at_pit[10] = 200.  # juice it up again
    >>> _raise_lake_to_limit(10, lake_map, area_map, master_pit_q,
    ...                      lake_q_dict, init_water_surface, lake_water_level,
    ...                      accum_area_at_pit, vol_rem_at_pit,
    ...                      accum_K_at_pit, lake_is_full, lake_spill_node,
    ...                      water_vol_balance_terms, neighbors, closednodes,
    ...                      drainingnodes, water_out_dict)

    Here, rise continues until the combined lake hits the drainingnode at 29.
    Note it doesn't feel the lower, closed node at 34:

    >>> np.all(np.equal(lake_map, np.array([-1, 43, 43, -1, 43, -1,
    ...                                     43,  7,  7, 43,  7, 43,
    ...                                     43,  7,  7,  7,  7, 43,
    ...                                     43,  7, 43,  7, 43, -1,
    ...                                     -1, 43, 43,  7,  7, 79,
    ...                                     -1, -1, -1, 43, -1, -1])))
    True
    >>> lake_q_dict[7].tasks_currently_in_queue()
    array([11,  6, 33,  4,  9, 26,  1,  2, 17, 22, 20, 12, 25, 18])
    >>> lake_water_level == {7: -1.0}
    True
    >>> accum_area_at_pit == {7: 44.0}
    True
    >>> vol_rem_at_pit == {7: 182.0}
    True
    >>> lake_is_full == {7: False}
    True
    >>> lake_spill_node == {7: 29}
    True
    """
    # We work upwards in elev from the current level, raising the level to
    # the next lowest node. We are looking for the first sign of flow
    # going into a node that isn't already in the queue that is *downhill*
    # of the current elevation. The idea here is that either a node is
    # accessible from below via another route, in which case it's already in
    # the queue, or it's over a saddle, and thus belongs to a separate pit.
    cpit = current_pit  # shorthand
    # First, check it's valid to call this in the first place. If this lake
    # already has a spill, then the lake is already raise, so *this* func
    # need not do anything. BUT we do need to route any excess onwards.
    if lake_spill_node[cpit] != -1:
        return

    tobreak = False  # delayed break flag
    nnodes = init_water_surface.size
    pit_nghb_code = cpit + nnodes
    # ^this code used to flag possible nodes to inundate, but that aren't
    # flooded yet.
    spill_code = pit_nghb_code + nnodes
    # ^more funky flagging to indicate a future spill
    surf_z = init_water_surface
    # We're calling this in the first place, so...
    lake_is_full[cpit] = False

    # This func assumes that when it receives a lake_q_dict, it will be pre-
    # loaded with (first, though this should be clear from the priority) the
    # cpit, or the node from which the iteration should start (i.e., the
    # node where the lake terminated last time), and then also that node's
    # immediate neighbors. This means this works:
    cnode = lake_q_dict[cpit].pop_task()

    while not tobreak:
        # print('New loop. lake_q', lake_q_dict[cpit].tasks_currently_in_queue())
        nnode = lake_q_dict[cpit].pop_task()
        if 0 <= lake_map[nnode] < nnodes:
            # this is not a valid next node, as it's already in a lake.
            # this can arise after lakes are merged...
            continue
        # print('Start at', cnode, 'Next', nnode, 'Code next', lake_map[nnode])
        # note that this is the next node we will try to flood, not the one
        # we are currently flooding. A node is incorporated into the lake at
        # the point when _raise_water_level is run.
        # cnode is the current node we are flooding.
        # is the cnode truly virgin?
        freshnode = (lake_map[cnode] < 0) or (
            lake_map[cnode] >= nnodes)
        # Now, try to raise the level to this height, without running out
        # of water...
        z_increment = surf_z[nnode] - lake_water_level[cpit]
        (filled, z_increment) = _raise_water_level(
            cpit, cnode, area_map, z_increment, lake_map,
            water_vol_balance_terms, accum_area_at_pit, vol_rem_at_pit,
            accum_K_at_pit, lake_spill_node,
            flood_from_zero=freshnode)
        lake_water_level[cpit] += z_increment
        if not filled:
            # print('Ran out of water:', cnode)
            # Now, the lake can't grow any more! (...at the moment.)
            # there's still potential here for this lake to grow more, and
            # if it does, the first node to rise will be this one. So,
            # stick the node back in the local lake queue, along with all
            # the higher neighbors that never got raised to before we severed
            # the iteration:
            lake_is_full[cpit] = True
            lake_q_dict[cpit].add_task(nnode, priority=surf_z[nnode])
            lake_q_dict[cpit].add_task(cnode, priority=lake_water_level[cpit])
            break

        # Now, a special case where we've filled the lake to this level,
        # but it turns out the next node is in fact the spill of an
        # adjacent lake...
        if lake_map[nnode] >= 2 * nnodes:
            # this condition only triggered (I think!) if we are now
            # raised to the spill level of an existing lake.
            # this then becomes break case (b) - except we cheat, and
            # don't actually break, since we're safe to continue raising
            # the new, composite lake as if it were a new one.
            cpit = _merge_two_lakes(
                this_pit=cpit, that_pit=lake_map[nnode] - 2*nnodes,
                z_topo=init_water_surface,
                lake_map=lake_map, lake_q_dict=lake_q_dict,
                lake_water_level=lake_water_level,
                accum_area_at_pit=accum_area_at_pit,
                vol_rem_at_pit=vol_rem_at_pit,
                accum_K_at_pit=accum_K_at_pit,
                lake_is_full=lake_is_full,
                lake_spill_node=lake_spill_node)
            # print('Merging!', cpit)
            cnode = nnode
            # print('Vol, A', vol_rem_at_pit[cpit], accum_area_at_pit[cpit])
            pit_nghb_code = cpit + nnodes
            spill_code = pit_nghb_code + nnodes
            continue
            # we allow the loop to continue after this, provided the
            # queues and dicts are all correctly merged
            # in this case, leave the cnode as it is... they're all at the
            # same level now
            # allow a nghb check just in case we have a super funky geometry,
            # but very likely to find new neighbors
        if drainingnodes[nnode]:
            lake_map[nnode] = spill_code
            lake_spill_node[cpit] = nnode
            try:
                water_out_dict[nnode] += vol_rem_at_pit[cpit]
                # this necessary in case there is already water spilling,
                # and more is added later from above (a cascade)
            except KeyError:
                water_out_dict[nnode] = vol_rem_at_pit[cpit]
            break
        cnode = nnode
        # ^Note, even in the merging case, we leave the actual sill as the next
        # node, so the current z_surf makes sense next time around

        # OK, so we've filled that node. Now let's consider where next...
        cnghbs = neighbors[cnode]
        # note that an "open" neighbor should potentially include the spill
        # node of an adjacent lake. This case dealt with in the raise step.
        for n in cnghbs:
            if not closednodes[n]:
                if lake_map[n] not in (cpit, pit_nghb_code):
                    # ^virgin node (n==-1), or unclaimed nghb of other lake
                    # (0<=n<nnodes and n!=cpit), or future spill of another
                    # lake are all permitted.
                    # Now, if we've got to here, then all open neighbors must
                    # lead upwards; a spur rather than a sill will have already
                    # been explored from below. So if the topo falls, current
                    # node is a true sill
                    if surf_z[n] < surf_z[cnode]:
                        # print('Down!', cnode, n, lake_water_level[cpit])
                        lake_map[cnode] = spill_code
                        # "magic" coding that lets this lake put "dibs" on that
                        # spill node.
                        lake_spill_node[cpit] = cnode
                        master_pit_q.add_task(
                            cpit, priority=lake_water_level[cpit])
                        tobreak = True
                        # ...but nevertheless, we want to keep loading the
                        # others so we can continue to use this list later if
                        # the lake acquires a new flux
                    else:
                        # print('Up!', cnode, n, lake_water_level[cpit])
                        lake_q_dict[cpit].add_task(n, priority=surf_z[n])
                        lake_map[n] = max((pit_nghb_code, lake_map[n]))
                        # note outlet codes take priority


def _merge_two_lakes(this_pit, that_pit, z_topo, lake_map,
                     lake_q_dict,
                     lake_water_level,
                     accum_area_at_pit,
                     vol_rem_at_pit,
                     accum_K_at_pit,
                     lake_is_full,
                     lake_spill_node):
    """
    Take two lakes that are known to be at the same level & in contact,
    and merge them.

    Remember, at the time of merge, any sill is NOT part of the merged lake
    (a call of _raise_water_level on the sill is required to actually
    inundate it).

    Parameters
    ----------
    this_pit : int
        The pit ID of the first lake.
    that_pit : int
        The pit ID of the second lake.
    z_topo : array
        The original topography of the surface on which the lakes occur, i.e.,
        without water. This is used to ensure the assigned lake ID corresponds
        to the deepest lake.
    lake_map : array of ints
        Grid-node array of -1s with IDs where a node is known to be flooded
        from that pit, and ID-related coding for lake adjacent or lake outlet
        nodes.
    lake_q_dict : dict of StablePriorityQueues
        A dict of queues holding nodes the node which begins the iteration
        (typically the pit node), along with the first set of nodes to
        explore. This is required to remove the unused pit ID from the keys.
    lake_water_level : dict
        The current (transient) water level at each lake. This is required
        to remove the unused pit ID from the keys.
    accum_area_at_pit : dict
        The total surface area of each lake, updated through time. This is
        required to remove the unused pit ID from the keys.
    vol_rem_at_pit : dict
        The "excess" water volume available at each lake, that will be used to
        continue to raise the water level until exhausted. This is required to
        remove the unused pit ID from the keys.
    accum_K_at_pit : dict
        Stores the current total K term in the (optional) water balance
        equation across the sum of all nodes in each lake. Needed to perform
        water balance as the lake level rises during the fill. This is
        required to remove the unused pit ID from the keys.
    lake_is_full : dict
        Stores the status of each distinct pit in the depression. This is
        required to remove the unused pit ID from the keys.
    lake_spill_node : dict
        the current spill node of each lake. -1 if not defined. This is
        required to remove the unused pit ID from the keys.

    Returns
    -------
    cpit : int
        ID of the lake as defined by the lower of the two pits provided.
    """
    # Note we *will* still have stuff in the master_pit_q that has the ID of
    # a replaced lake. Deal with this by spotting that those lake codes no
    # longer work as dict keys...
    # Check they've got to the same level (really ought to have!)
    assert np.isclose(lake_water_level[this_pit], lake_water_level[that_pit])
    nnodes = lake_map.size
    if z_topo[this_pit] < z_topo[that_pit]:
        low_pit = this_pit
        hi_pit = that_pit
    else:
        low_pit = that_pit
        hi_pit = this_pit

    # merge the queues
    # print('this queue', lake_q_dict[this_pit].tasks_currently_in_queue())
    # print('that queue', lake_q_dict[that_pit].tasks_currently_in_queue())
    lake_q_dict[low_pit].merge_queues(lake_q_dict[hi_pit])
    # print('new queue', lake_q_dict[low_pit].tasks_currently_in_queue())
    _ = lake_q_dict.pop(hi_pit)
    # merge the maps
    lake_map[lake_map == hi_pit] = low_pit
    lake_map[lake_map == (hi_pit + nnodes)] = low_pit + nnodes
    # merge the various dicts
    for dct in (accum_area_at_pit, vol_rem_at_pit, accum_K_at_pit):
        dct[low_pit] += dct.pop(hi_pit)
    for dct in (lake_water_level, lake_is_full, lake_spill_node):
        _ = dct.pop(hi_pit)
    lake_spill_node[low_pit] = -1  # just in case

    return low_pit


def _get_float_water_vol_balance_terms(water_vol_balance_terms,
                                       area_map, node):
    """
    Takes the floats or arrays of the water_vol_balance_terms, and returns
    the appropriate float value at the given node, already having termed a
    "loss per unit depth" into a "volume loss".

    Examples
    --------
    >>> import numpy as np
    >>> area_map = 2. * np.arange(4) + 1.
    >>> c, K = _get_float_water_vol_balance_terms((-3., np.arange(4)),
    ...                                           area_map, 0)
    >>> type(c) is np.float64
    True
    >>> type(K) is np.float64
    True
    >>> np.isclose(c, -3.)
    True
    >>> np.isclose(K, 0.)
    True

    >>> c, K = _get_float_water_vol_balance_terms((-1 * np.arange(4), 2.),
    ...                                           area_map, 3)
    >>> type(c) is np.float64
    True
    >>> type(K) is np.float64
    True
    >>> np.isclose(c, -21.)
    True
    >>> np.isclose(K, 14.)
    True
    """
    cell_A = area_map[node]
    if type(water_vol_balance_terms[0]) is np.ndarray:
        c = water_vol_balance_terms[0][node]
    else:
        c = water_vol_balance_terms[0]
    if type(water_vol_balance_terms[1]) is np.ndarray:
        K = water_vol_balance_terms[1][node]
    else:
        K = water_vol_balance_terms[1]
    return (c*cell_A, K*cell_A)


def _raise_water_level(cpit, cnode, area_map, z_increment,
                       lake_map, water_vol_balance_terms,
                       accum_area_at_pit, vol_rem_at_pit, accum_K_at_pit,
                       lake_spill_node, flood_from_zero=True):
    """
    Lift water level from the foot of a newly flooded node surface to the
    level of the next lowest, while honouring and water balance losses.

    Note: does not actually raise the water_level value! Returns the
    change in lake elevation instead.

    Parameters
    ----------
    cpit : int
        Current pit identifying this lake.
    cnode : int
        The node that is currently being inundated.
    area_map : array of floats
        Grid-node array of the areas of each cell; typically the output of
        a grid cell_area_at_node property.
    z_increment : float
        The total change in elevation between the current level and the target
        level.
    lake_map : array of ints
        Grid-node array of -1s with IDs where a node is known to be flooded
        from that pit, and ID-related coding for lake adjacent or lake outlet
        nodes.
    water_vol_balance_terms : (float c, float K)
        polyparams for a linear fit of the water vol balance func, such that
        dV = A * (K*z + c). Fit must be linear.
    accum_area_at_pit : dict
        The total surface area of each lake, updated through time.
    vol_rem_at_pit : dict
        The "excess" water volume available at each lake, that will be used to
        continue to raise the water level until exhausted.
    accum_K_at_pit : dict
        Stores the current total K term in the (optional) water balance
        equation across the sum of all nodes in each lake. Needed to perform
        water balance as the lake level rises during the fill.
    lake_spill_node : dict
        the current spill node of each lake. -1 if not defined.
    flood_from_zero : bool
        Flag to indicate whether the fill was from a node already inundated,
        or the inundation of a fresh node. (established from lake_map codes)

    Returns
    -------
    (full_fill, z_increment) : (bool, float)
        Did the node fully fill, and what increment of water depth was added
        in the end?

    Examples
    --------

    >>> import numpy as np
    >>> lake_map_init = np.ones(9, dtype=int) * -1
    >>> lake_map_init[4] = 5
    >>> lake_map = lake_map_init.copy()
    >>> area_map = np.array([ 2., 2., 2.,
    ...                       4., 3., 2.,
    ...                       1., 2., 2.])
    >>> area_map[3] = 4.
    >>> area_map[4] = 3.
    >>> accum_area_at_pit = {5: 3.}
    >>> vol_rem_at_pit = {5: 100.}
    >>> accum_K_at_pit = {5: 0.}
    >>> lake_spill_node = {5: 8}
    >>> full, dz = _raise_water_level(cpit=5, cnode=3, area_map=area_map,
    ...                               z_increment=3., lake_map=lake_map,
    ...                               water_vol_balance_terms=(0., 0.),
    ...                               accum_area_at_pit=accum_area_at_pit,
    ...                               vol_rem_at_pit=vol_rem_at_pit,
    ...                               accum_K_at_pit=accum_K_at_pit,
    ...                               lake_spill_node=lake_spill_node,
    ...                               flood_from_zero=True)
    >>> full
    True
    >>> dz == 3.
    True
    >>> np.all(np.equal(lake_map, np.array([-1, -1, -1,
    ...                                       5,  5, -1,
    ...                                      -1, -1, -1])))
    True
    >>> np.isclose(accum_area_at_pit[5], 7.)
    True
    >>> np.isclose(vol_rem_at_pit[5], 79.)
    True
    >>> np.isclose(accum_K_at_pit[5], 0.)
    True
    >>> lake_spill_node[5]
    8

    >>> full, dz = _raise_water_level(cpit=5, cnode=5, area_map=area_map,
    ...                               z_increment=4., lake_map=lake_map,
    ...                               water_vol_balance_terms=(-5., 0.),
    ...                               accum_area_at_pit=accum_area_at_pit,
    ...                               vol_rem_at_pit=vol_rem_at_pit,
    ...                               accum_K_at_pit=accum_K_at_pit,
    ...                               lake_spill_node=lake_spill_node,
    ...                               flood_from_zero=True)
    >>> full
    True
    >>> dz == 4.
    True
    >>> np.all(np.equal(lake_map, np.array([-1, -1, -1,
    ...                                      5,  5,  5,
    ...                                     -1, -1, -1])))
    True
    >>> np.isclose(accum_area_at_pit[5], 9.)
    True
    >>> np.isclose(vol_rem_at_pit[5], 33.)
    True
    >>> np.isclose(accum_K_at_pit[5], 0.)
    True
    >>> lake_spill_node[5]
    8

    >>> full, dz = _raise_water_level(cpit=5, cnode=6, area_map=area_map,
    ...                               z_increment=1., lake_map=lake_map,
    ...                               water_vol_balance_terms=(-34., -3.),
    ...                               accum_area_at_pit=accum_area_at_pit,
    ...                               vol_rem_at_pit=vol_rem_at_pit,
    ...                               accum_K_at_pit=accum_K_at_pit,
    ...                               lake_spill_node=lake_spill_node,
    ...                               flood_from_zero=True)
    >>> full
    False
    >>> dz == 0.
    True
    >>> np.all(np.equal(lake_map, np.array([-1, -1, -1,
    ...                                      5,  5,  5,
    ...                                      5, -1, -1])))
    True
    >>> np.isclose(accum_area_at_pit[5], 10.)
    True
    >>> np.isclose(vol_rem_at_pit[5], -1.)
    True
    >>> np.isclose(accum_K_at_pit[5], 0.)
    True
    >>> lake_spill_node[5]
    -1
    >>> lake_spill_node[5] = 8
    >>> vol_rem_at_pit[5] = 87.

    >>> full, dz = _raise_water_level(cpit=5, cnode=6, area_map=area_map,
    ...                               z_increment=2., lake_map=lake_map,
    ...                               water_vol_balance_terms=(-34., -3.),
    ...                               accum_area_at_pit=accum_area_at_pit,
    ...                               vol_rem_at_pit=vol_rem_at_pit,
    ...                               accum_K_at_pit=accum_K_at_pit,
    ...                               lake_spill_node=lake_spill_node,
    ...                               flood_from_zero=False)
    >>> full
    True
    >>> dz == 2.
    True
    >>> np.all(np.equal(lake_map, np.array([-1, -1, -1,
    ...                                      5,  5,  5,
    ...                                      5, -1, -1])))
    True
    >>> np.isclose(accum_area_at_pit[5], 10.)  # ...still
    True
    >>> np.isclose(vol_rem_at_pit[5], 61.)  # -6 loss, -20 rise
    True
    >>> np.isclose(accum_K_at_pit[5], -3.)
    True
    >>> lake_spill_node[5]
    8

    >>> full, dz = _raise_water_level(cpit=5, cnode=2, area_map=area_map,
    ...                               z_increment=10., lake_map=lake_map,
    ...                               water_vol_balance_terms=(
    ...                                 np.arange(9), np.ones(9)),
    ...                               accum_area_at_pit=accum_area_at_pit,
    ...                               vol_rem_at_pit=vol_rem_at_pit,
    ...                               accum_K_at_pit=accum_K_at_pit,
    ...                               lake_spill_node=lake_spill_node,
    ...                               flood_from_zero=True)
    >>> full
    False
    >>> np.isclose(dz, 5.)
    True
    >>> np.all(np.equal(lake_map, np.array([-1, -1,  5,
    ...                                      5,  5,  5,
    ...                                      5, -1, -1])))
    True
    >>> np.isclose(accum_area_at_pit[5], 12.)
    True
    >>> np.isclose(vol_rem_at_pit[5], 0.)  # +4, then -5 (=-3+2) loss, -60 rise
    True
    >>> np.isclose(accum_K_at_pit[5], -3.)  # node not filled, so not changed
    True
    >>> lake_spill_node[5]
    -1
    """
    # first up, one way or another this is now in the lake:
    lake_map[cnode] = cpit
    if flood_from_zero:
        accum_area_at_pit[cpit] += area_map[cnode]
    V_increment_to_fill = (
        z_increment * accum_area_at_pit[cpit])
    c_added_at_start, K_to_lift = _get_float_water_vol_balance_terms(
        water_vol_balance_terms, area_map, cnode)

    if flood_from_zero:
        vol_rem_at_pit[cpit] += c_added_at_start
        if vol_rem_at_pit[cpit] < 0.:
            lake_spill_node[cpit] = -1
            return (False, 0.)  # return a zero, but the node is still flooded

    V_gain_in_full_fill = (
        K_to_lift + accum_K_at_pit[cpit]) * z_increment  # likely <= 0.
    V_increment_to_fill -= V_gain_in_full_fill  # bigger if lossy
    # note that it's very possible that we fill the node then just add to
    # the total if the pit is, in fact, gaining... (careful not to double
    # account any gains wrt the flow routing)
    if vol_rem_at_pit[cpit] < V_increment_to_fill:
        # we don't have the water to get onto the next level
        frac = vol_rem_at_pit[cpit]/V_increment_to_fill
        lake_spill_node[cpit] = -1
        z_available = frac * z_increment
        vol_rem_at_pit[cpit] = 0.
        return (False, z_available)
    else:
        # successful fill of node
        # increment the K term, so when we fill next this node is already
        # reflected:
        accum_K_at_pit[cpit] += K_to_lift
        vol_rem_at_pit[cpit] -= V_increment_to_fill
        return (True, z_increment)


def _route_outlet_to_next(outlet_ID, lake_at_head, flow__receiver_nodes,
                          link_lengths, links_at_node, z_surf, lake_map,
                          neighbors, closednodes, drainingnodes):
    """
    Take an outlet node, and then follow the steepest descent path until it
    reaches a distinct lake.  Transmission losses are here NOT accounted for;
    we assume it all comes out in the wash given they likely were applied
    during the initial flow accumulation that provided the inputs to the pits.

    Note that this function does not actually route the flow; it simply
    finds and identifies the next lake or outlet downstream. Do this
    discharge rewiring elsewhere.

    Parameters
    ----------
    outlet_ID : int
        The outlet ID of the lake from which flow is being routed.
    lake_at_head : int
        The lake ID of the lake from which flow is being routed.
    flow__receiver_nodes : array of ints
        The nodes already known to be the steepest descent path on the
        unflooded terrain. Typically the output from a many-to-one
        FlowDirector.
    link_lengths : array of floats
        The length of each link on the grid. Typically grid.length_of_link.
    links_at_node : (nnodes, max_nlinks_at_node)
        The links at each node on the grid, -1 where undefined. Typically
        grid.links_at_node or d8 equivalent.
    z_surf : array of floats
        The topography on the grid. Only operated upon where the terrain is
        not flooded, so this can be either the flooded or unflooded data.
    lake_map : array of ints
        Grid-node array of -1s with IDs where a node is known to be flooded
        from that pit, and ID-related coding for lake adjacent or lake outlet
        nodes.
    neighbors : (nnodes, max_nneighbors) array
        The neighbors at each node on the grid. Typically obtained direct
        from a grid, as the adjacent_nodes_at_node or the equivalent
        d8 structure.
    closednodes : array of bool
        Nodes not to be explored by the algorithm.
    drainingnodes : array of bool
        Once a lake raises itself to the level of one of these nodes, it will
        cease to rise any more, ever (i.e., these are boundary nodes!)

    Returns
    -------
    (next_lake_ID, draining_node_out) : (int, int)
        next_lake_ID is the ID of the next lake reached by flow from this
        spill point. -1 if flow drains into a drainingnode. draining_node_out
        is the node ID of one of the drainingnode which the flow reaches, -1
        if the flow terminates in a lake. i.e., one of these values will always
        be -1.

    Examples
    --------
    >>> import numpy as np
    >>> from landlab import RasterModelGrid
    >>> z_surf = np.array([ 10.,  1., 10., 10., 10.,
    ...                      1.,  3.,  2.,  1.,  0.,
    ...                     10., 10.,  1., 10., 10.,])
    >>> lake_map = np.array([ 19, -1, -1, -1, 21,
    ...                        4, 34, -1, 21,  6,
    ...                       19, -1, -1, -1, 21])
    >>> spill = 6
    >>> lake_at_head = 4
    >>> closednodes = np.zeros(15, dtype=bool)
    >>> closednodes[[1, 12]] = True
    >>> drainingnodes = np.zeros(15, dtype=bool)
    >>> drainingnodes[12] = True
    >>> mg = RasterModelGrid((3, 5), (1., 2.))
    >>> flow_rec_node = np.array([ 5,  1,  7,  8,  9,
    ...                            5,  1, 12,  9,  9,
    ...                            5,  6, 12,  8,  9])
    >>> (lake, drain) = _route_outlet_to_next(
    ...     spill, lake_at_head, flow_rec_node, mg.length_of_link,
    ...     mg.links_at_node, z_surf, lake_map, mg.adjacent_nodes_at_node,
    ...     closednodes, drainingnodes)
    >>> lake
    6
    >>> drain
    -1
    >>> closednodes[12] = False
    >>> (lake, drain) = _route_outlet_to_next(
    ...     spill, lake_at_head, flow_rec_node, mg.length_of_link,
    ...     mg.links_at_node, z_surf, lake_map, mg.adjacent_nodes_at_node,
    ...     closednodes, drainingnodes)
    >>> lake
    -1
    >>> drain
    12
    """
    # first, resolve possible issue with FD of the outlet. It's possible it
    # currently drains back into the lake, which we need to correct.
    nghbs = neighbors[outlet_ID]
    nnodes = lake_map.size
    steepest_nghb = _steepest_valid_neighbor(
        outlet_ID, nghbs, lake_at_head, link_lengths, links_at_node, z_surf,
        lake_map, closednodes, drainingnodes)
    lake_code_steepest = lake_map[steepest_nghb]
    while (lake_code_steepest < 0 or lake_code_steepest >= nnodes):
        if drainingnodes[steepest_nghb]:
            return (-1, steepest_nghb)
        new_steepest_nghb = flow__receiver_nodes[steepest_nghb]
        if closednodes[new_steepest_nghb]:
            # locally forbidden, not necessarily seen by flow routing
            # so redirect here
            new_steepest_nghb = _steepest_valid_neighbor(
                steepest_nghb, neighbors[steepest_nghb], lake_at_head,
                link_lengths, links_at_node, z_surf, lake_map, closednodes,
                drainingnodes)
        steepest_nghb = new_steepest_nghb
        lake_code_steepest = lake_map[steepest_nghb]
    return (lake_code_steepest, -1)


def _steepest_valid_neighbor(node, nghbs, lake_at_head, link_lengths,
                             links_at_node, z_surf, lake_map, closednodes,
                             drainingnodes):
    """
    Find the node at the end of the steepest downward link from a node, that
    is not a closednode, nor in the lake at the head of the flow.

    Parameters
    ----------
    node : int
        The node in the flow path at which to determine the steepest downslope
        neighbor.
    nghbs : array of ints
        The neighbors at the node. -1s where undefined are fine.
    lake_at_head : int
        The lake ID of the lake at the head of the flow path.
    link_lengths : array of floats
        The length of each link on the grid. Typically grid.length_of_link.
    links_at_node : (nnodes, max_nlinks_at_node)
        The links at each node on the grid, -1 where undefined. Typically
        grid.links_at_node or d8 equivalent.
    z_surf : array of floats
        The topography on the grid. Only operated upon where the terrain is
        not flooded, so this can be either the flooded or unflooded data.
    lake_map : array of ints
        Grid-node array of -1s with IDs where a node is known to be flooded
        from that pit, and ID-related coding for lake adjacent or lake outlet
        nodes.
    closednodes : array of bool
        Nodes not to be explored by the algorithm.
    drainingnodes : array of bool
        Once a lake raises itself to the level of one of these nodes, it will
        cease to rise any more, ever (i.e., these are boundary nodes!)

    Returns
    -------
    steepest_nghb : int
        The node ID of the valid adjacent node with the steepest downward
        slope from node.

    Examples
    --------

    >>> import numpy as np
    >>> from landlab import RasterModelGrid
    >>> z_surf = np.array([ 10.,  1., 10., 10., 10.,
    ...                      1.,  3.,  2.,  1.,  0.,
    ...                     10., 10.,  1., 10., 10.,])
    >>> lake_map = np.array([ 19, -1, -1, -1, 21,
    ...                        4, 34, -1, 21,  6,
    ...                       19, -1, -1, -1, 21])
    >>> closednodes = np.zeros(15, dtype=bool)
    >>> closednodes[1] = True
    >>> drainingnodes = np.zeros(15, dtype=bool)
    >>> drainingnodes[12] = True
    >>> mg = RasterModelGrid((3, 5), (1., 2.))
    >>> for node in (6, 7, 12, 14):
    ...     nghbs = mg.adjacent_nodes_at_node[node]
    ...     print(_steepest_valid_neighbor(node, nghbs, 4, mg.length_of_link,
    ...                                    mg.links_at_node, z_surf, lake_map,
    ...                                    closednodes, drainingnodes))
    7
    12
    12
    9
    """
    if drainingnodes[node]:
        return node  # this really should be redundant in a sensible use case
    nnodes = lake_map.size
    lake_code_nghbs = lake_map[nghbs]
    valid_positions = np.logical_and.reduce((
        np.not_equal(nghbs, -1), np.not_equal(closednodes[nghbs], True),
        np.not_equal(lake_code_nghbs, lake_at_head)))
    valid_nghbs = nghbs[valid_positions]
    slopes = ((z_surf[valid_nghbs] - z_surf[node]) /
              link_lengths[links_at_node[node, valid_positions]])
    steepest_nghb = valid_nghbs[np.argmin(slopes)]
    return steepest_nghb


class LakeFillerWithFlux(LakeMapperBarnes):
    """
    """

    def run_one_step(dt):
        """
        """
        # First, we need a conventional lake fill. This involves calling the
        # FillSinksBarnes component. Nice, FAST special cases where lakes fill
        # up or only have one pit, which we can use to accelerate things,
        # so this is very much worth it.
        pass

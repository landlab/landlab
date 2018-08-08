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


class LakeEvaporator(Component):
    """
    Component to take a water surface overlying a bedrock surface as a series
    of (specified) existing lakes, then reduce the volume of water in each
    lake (i.e., drop the water surface) according to a loss term specified
    at each lake node.

    The component operates by knowing the maximum fill condition that a lake
    could reach, then


    Parameters
    ----------
    grid : ModelGrid
        The grid on which to operate the component.
    water_surface : array or field name
        The free surface of the topography.
    ground_surface : array or field name
        The (rock, soil) surface at or underlying the water_surface.
    water_balance : array or field name
        The driving water balance term, negative for loss, positive for gain.
        Expressed as flux in each cell, i.e., [L/T]
        Note that any gain that would lift a lake beyond full is considered
        lost; it is NOT transferred to the next lake downstream.
    lake_dict : dict of {outlet: (iterable of lake nodes)}
        A lake dictionary providing a description of the lakes across the
        landscape **if every pit in the landscape were to be completely
        filled**. This is required for efficiency of the component.
    """
    def __init__(self, grid, water_surface, ground_surface, water_balance,
                 lake_dict):
        """
        Initialise the component.
        """
        self._grid = grid
        self._water = return_array_at_node(water_surface)
        self._ground = return_array_at_node(ground_surface)
        assert np.all(np.less_equal(self._ground, self._water)), \
            "Ground surface must be at or below the water surface everywhere!"
        self._water_balance = return_array_at_node(water_balance)
        self._lakes = lake_dict

    def _solve_one_lake(self, lake_out, lake_node_array):
        """
        """
        total_water_balance_vol = (
            self._water_balance[lake_node_array] *
            self.grid.cell_area_at_node[lake_node_array]).sum()
        if np.isclose(total_water_balance_vol, 0.):
            pass  # no adjustment here
        else:
            init_water_level = self._water[lake_node_array[0]]
            max_level = self._ground[lake_out]
            # check the fill level makes sense... Also implicit in the method
            assert np.allclose(self._water[lake_node_array], init_water_level)
            # sort the nodes in the lake by ground elev
            lake_node_elevs = self._ground[lake_node_array]
            node_order = np.argsort(lake_node_elevs)
            elevs_in_order = self._ground[node_order]
            nodes_in_order = lake_node_array[node_order]
            areas_in_order = self.grid.cell_area_at_node[nodes_in_order]
            elev_diffs_in_order = np.ediff1d(elevs_in_order)
            accumulated_area_in_order = np.cumsum(areas_in_order)
            vol_diffs_in_order = (
                elev_diffs_in_order * accumulated_area_in_order[:-1])
            # ^this is the volume required to change the volume by one
            # "elevation step" at the various heights
            total_vol_in_order = np.cumsum(vol_diffs_in_order)
            # so, where does the starting elev appear in the node order?
            add_index = np.searchsorted(elevs_in_order, init_water_level)
            # Note we *don't* need to know the init vol. Doing it relative
            # to first position should work...
            if add_index - 2 < 0:
                first_vol = 0.
                assert total_water_balance_vol > 0.
                slc = slice(add_index, None, 1)  # we will advance up
            else:
                vol_to_above = np.amin((total_vol_in_order[add_index-1], 0.))
                vol_to_below = np.amin((total_vol_in_order[add_index-2], 0.))
                # if bone dry pit, add_index==0, both are 0.
                elev_above = elevs_in_order[add_index]
                elev_below = elevs_in_order[add_index - 1]

                # now whether we increase or drop matters:
                vol_first_slab = vol_to_above - vol_to_below
                if total_water_balance_vol > 0.:
                    slc = slice(add_index, None, 1)  # we will advance up
                    elev_diff_above = elevs_in_order[add_index] - init_water_level
                    elev_diff = elev_diff_above
                else:
                    slc = slice(add_index-1, None, -1)  # work down
                    elev_diff_below = init_water_level - elevs_in_order[add_index]
                    elev_diff = elev_diff_below
                if
                first_vol = vol_first_slab * elev_diff/(elev_above-elev_below)



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
    lake_water_volume_balance = {}  # tracks any total water balance on lake
    accum_Ks_at_pit = {}
    lake_is_full = {}
    for pit, vol, A, level in zip(pit_nodes_in_order, Vw_at_pits_in_order,
                                  area_map[pit_node_sort],
                                  zwater_nodes_in_order):
        vol_rem_at_pit[pit] = vol
        accum_area_at_pit[pit] = A
        lake_water_level[pit] = level
        lake_is_full[pit] = False
        lake_water_volume_balance[pit] = 0.
        accum_Ks_at_pit[pit] = 0.
    # these will get adjusted thru fill process
    # ^ idea here w vol_rem is that we will propagate the info on what the
    # total V assoc w each pit as well as the pit ID itself. This will let us
    # adjust the total vol balance as the lake spreads (e.g., if infiltration
    # matters)

    # note there's a special case hiding in here where one subbasin can fill
    # to a sill, but the fill in the next depression isn't enough to reach
    # it. Clearly we need to transfer this overtop into the next basin, but
    # not clear how to establish this connectivity...

    flood_order = StablePriorityQueue()
    for pit_ID in pit_nodes_in_order:
        flood_order.add_task(pit, priority=water_z[pit])

    tomerge_flag = -1  # if this has a +ve number, will merge lakes
    lakes_weve_look_at = set()
    lakes_to_work = set()
    # distinguish these so we can't put lakes we've filled back in the list
    while True:
        # the trick here is to NOT actually update the elevs at the nodes.
        # instead we simply move the "virtual" water level up, so we don't
        # have to explicitly track the levels at all the nodes, and modify
        # them every step
        # grab the next node:
        current_node = flood_order.pop_task()
        current_node_pit = lake_map[current_node]
        if maybe_pit not in lakes_weve_look_at:
            lakes_weve_look_at.append(current_node)
            lakes_to_work.append(current_node)
        node_neighbors = neighbors_at_nodes[current_node]
        current_node_zsurf = lake_water_level[current_node_pit]

        # load up the queue. Likely one of these will be next, so do it 1st
        for nghb in np.nditer(node_neighbors):
            if lake_map[nghb] < 0 or lake_map[nghb] += num_pits:  # virgin territory
                # nodes on the perim of other pits are fair game!
                zsurf_nghb = water_z[nghb]  # actual water level at start

#### START HERE






                if zsurf_nghb >= current_node_zsurf:
                    if zsurf_nghb > nextlowest_zsurf:
                        flood_order.add_task(nghb, priority=zsurf_nghb)
                        lake_map[nghb] = currentnode_pit + num_pits
                        # this is clever flagging to distinguish which lake
                        # a perimeter node is part of, but that it isn't
                        # yet flooded...
                    else:
                        flood_order.add_task(nextlowest,
                                             priority=nextlowest_zsurf)
                        nextlowest = nghb
                        nextlowest_zsurf = zsurf_nghb
            else:
                ######FOUND A SILL, TREAT HERE...? -> YES
                # what is the lake we found?
                found_lake = lake_map[nghb]
                # flag it for merge, but we can only
                # Think this will be easier if we merge this lake into the
                # one we're working on right now...
                lake_map[lake_map == found_lake] = currentnode_pit
                # merge the info in the dicts...
                for dict_to_merge in (vol_rem_at_pit, accum_area_at_pit,)
                    dict_to_merge[currentnode_pit] += dict_to_merge.pop(
                        found_lake)

                lake_water_volume_balance
                accum_Ks_at_pit
                lake_is_full


                ## lake_water_level[currentnode_pit] = lake_water_level[found_lake]
                # not yet. Still need to raise the currentnode to the level.
                # or not at all? Still need to budget all this raising...





                pass

        # we're going to need to work lakes simulateously, so:
        for currentnode_pit in lakes_to_work:      ########CALL THIS PIT
            if lake_is_full[currentnode_pit]:
                # stop treating this pit entirely
                lakes_to_work.remove(currentnode_pit)
                #### do we need to take out the perim nodes that might remain?
                continue
            # for the purposes of the logic, a flooded node (i.e., tagged w ID)
            # is one which is contiguous w a filled node. Does NOT have to be
            # filled to its level
            zsurf_node = lake_water_level[currentnode_pit]  # formerly water_z
            # at this point we know the water is lapping at the foot of this node
            node_neighbors = neighbors_at_nodes[current_node]
            currentnode_A = accum_area_at_pit[currentnode_pit]
            currentnode_V_rem = vol_rem_at_pit[currentnode_pit]
            # both already incremented at end of last step

            # what's the next node likely to be right now?
            nextlowest = flood_order.pop_task()
            nextlowest_zsurf = lake_water_level[nextlowest]

            # load up the queue. Likely one of these will be next, so do it 1st
            for nghb in np.nditer(node_neighbors):
                if lake_map[nghb] < 0 or lake_map[nghb] += num_pits:  # virgin territory
                    # nodes on the perim of other pits are fair game!
                    zsurf_nghb = water_z[nghb]
                    if zsurf_nghb >= zsurf_node:
                        if zsurf_nghb > nextlowest_zsurf:
                            flood_order.add_task(nghb, priority=zsurf_nghb)
                            lake_map[nghb] = currentnode_pit + num_pits
                            # this is clever flagging to distinguish which lake
                            # a perimeter node is part of, but that it isn't
                            # yet flooded...
                        else:
                            flood_order.add_task(nextlowest,
                                                 priority=nextlowest_zsurf)
                            nextlowest = nghb
                            nextlowest_zsurf = zsurf_nghb
                else:
                    ######FOUND A SILL, TREAT HERE...? -> YES
                    # what is the lake we found?
                    found_lake = lake_map[nghb]
                    # flag it for merge, but we can only
                    # Think this will be easier if we merge this lake into the
                    # one we're working on right now...
                    lake_map[lake_map == found_lake] = currentnode_pit
                    # merge the info in the dicts...
                    for dict_to_merge in (vol_rem_at_pit, accum_area_at_pit,)
                        dict_to_merge[currentnode_pit] += dict_to_merge.pop(
                            found_lake)

                    lake_water_volume_balance
                    accum_Ks_at_pit
                    lake_is_full


                    ## lake_water_level[currentnode_pit] = lake_water_level[found_lake]
                    # not yet. Still need to raise the currentnode to the level.
                    # or not at all? Still need to budget all this raising...





                    pass

            if lake_map[nextlowest] >= 0:
                # found an existing flooded surface, i.e., A SILL
                # Note possible for this to be a sill now when it wasn't before...
                # COME BACK AND SAVE THIS INFO
                # flag the sill, then route the discharge down into the next
                # depression?? (grab the flowdirs?)
                # do we mean continue? maybe still need to flood up. Or do later?
                # THINK THIS IS REDUNDANT SO
                raise AssertionError

            zsurf_next = water_z[nextlowest]
            z_increment = zsurf_next - zsurf_node
            # try to flood up to next
            V_increment_to_fill = (
                z_increment * currentnode_A -
                lake_water_volume_balance[currentnode_pit])
            # & remember to update the vol_balance below
            if V_increment_to_fill < 0.:
                # this really should not happen!!
                raise AssertionError, "lake is gaining too much direct input!!"
            # Note that the addition of the vol balance makes it harder to
            # actually fill the step in most cases.

            if currentnode_V_rem <= V_increment_to_fill:
                # we don't have the water to get onto the next level
                lake_is_full[currentnode_pit] = True
                z_available = currentnode_V_rem/currentnode_A
                lake_water_level[currentnode_pit] += z_available
                # no need to update water balance terms, since (a) this lake is
                # done, and (b) k's anc c's don't change anyway  as lake area is
                # static
                # Now, the lake can't grow any more! (...at the moment.)
            else:
                # successful fill of node
                vol_rem_at_pit[currentnode_pit] -= V_increment_to_fill
                lake_water_level[currentnode_pit] += z_increment
                # propagate the lake info into this next node, if it's
                # contiguous...
                if lake_map[current_node] == lake_map[nextlowest]:
                    # update the dicts describing the lake, as it floods onto
                    # the new surface.
                    _propagate_lake_to_next_node(currentnode_pit, nextlowest,
                                                 area_map, z_increment,
                                                 accum_area_at_pit,
                                                 lake_water_volume_balance,
                                                 water_vol_balance_terms,
                                                 accum_Ks_at_pit, c_to_add)
                else:
                    # NOT contiguous. In this case, we put it back in the
                    # queue, as it has more growing to do to reach an actual
                    # neighbor:
                    flood_order.add_task(
                        current_node, priority=lake_water_level[currentnode_pit])
                    # TODO: What is the surface to use at various pts? Check.
                    # note that due to the ordered queue, we will handle the
                    # two pits that must be dealt with here in an alternating
                    # fashion. This is important, since we are bringing both
                    # up to the same level simulataneously.
                    # Also, update the loss due to the increased depth:
                    lake_water_volume_balance[currentnode_pit] += (
                        accum_Ks_at_pit[currentnode_pit] * z_increment)
            # now, simply move to this next lowest node
            current_node = nextlowest


def _raise_lake_to_limit(grid, current_pit, lake_dict, lake_map, master_pit_q,
                         lake_q_dict, init_water_surface,
                         lake_water_level,
                         accum_area_at_pit,
                         vol_rem_at_pit,
                         lake_water_volume_balance,
                         accum_Ks_at_pit,
                         lake_is_full,
                         lake_spill_node,
                         water_vol_balance_terms,
                         neighbors, closednodes,
                         ):
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
    neighbors : (nnodes, max_nneighbors) array
        The neighbors at each node on the grid.
    closednodes : nodes not to be explored by the algorithm.
    Current pit : int
        Node ID of the pit that uniquely identifies the current lake that
        we are filling.
    lake_map : array of ints
        Grid-node array of -1s with IDs where a node is known to be flooded
        from that pit.
    master_pit_q: StablePriorityQueue
        This queue holds the current lakes in the grid, ordered by water
        surface elevation priority. The current lake will have just been
        popped off this list, exposing the next highest surface.
    lake_q_dict : dict of StablePriorityQueues
        A dict of queues holding nodes known to be on the perimeter of each
        lake. We hold this information since we return to each lake repeatedly
        to update the elevs. If the lake hasn't been filled before, this will
        just contain the pit node.
    init_water_surface : array
        The original topo of the water surface at the start of the step,
        i.e., before we start moving the lake water levels this time around.
    lake_water_level : dict
        The current (transient) water level at each lake.
     accum_area_at_pit : dict
        The total surface area of each lake
     vol_rem_at_pit : dict
        The "excess" water volume available at each lake, that will be used to
        continue to raise the water level until exhausted.
     lake_water_volume_balance : dict
        The volume loss expected from each lake under current lake area and
        depth.
     accum_Ks_at_pit : dict
        Stores the current total k term in the (optional) water balance
        equation across the sum of all nodes in each lake.
     lake_is_full : dict
        Stores the status of each distinct pit in the depression.
####probably redundant

     lake_spill_node : dict
        the current spill node of each lake. -1 if not defined.
    water_vol_balance_terms : (float c, float K)
        polyparams for a linear fit of the water vol balance func, such that
        V = A * (K*z + c)
    """
    # We work upwards in elev from the current level, raising the level to
    # the next lowest node. We are looking for the first sign of flow
    # going into a node that isn't already in the queue that is *downhill*
    # of the current elevation. The idea here is that either a node is
    # accessible from below via another route, in which case it's already in
    # the queue, or it's over a saddle, and thus belongs to a separate pit.
    tobreak = False  # delayed break flag
    cpit = current_pit  # shorthand
    pit_nghb_code = cpit + grid.number_of_nodes
    # ^this code used to flag possible nodes to inundate, but that aren't
    # flooded yet.
    spill_code = pit_nghb_code + grid.number_of_nodes
    # ^more funky flagging to indicate a future spill
    surf_z = init_water_surface

    ###### define cnode init here, and define neighbor queue correctly
    # might have to change ordering in the loop.
    cnode = cpit
    # add the neighbors here so they get immediately explored?
    # No, this probably needs doing **in the init steps (of the calling method)**
    # so everything is consistent
    # ...for that init of the method...
    # cnghbs = neighbors[cpit]
    # pit_nghb_code = cpit + grid.number_of_nodes
    # for n in np.nditer(cnghbs):
    #     if not closednodes[n]:
    #         # by definition, all neighbors of a pit are uphill, so
    #         lake_q_dict[cpit].add_task(n, priority=surf_z[n])
    #         lake_map[n] = pit_nghb_code
    #         # note that these potential neighbors can happily appear in more than
    #         # one list if they're spills... and indeed, this is implicit in how the
    #         # algorithm works!

    ####### we should optionally be saving all the volume information...?

    while not tobreak:
        nnode = lake_q_dict[cpit].pop_task()
        # note that this is the next node we will try to flood, not the one
        # we are currently flooding. A node is incorporated into the lake at
        # the point when _raise_water_level is run.
        # cnode is the current node we are flooding.
        # is the cnode truly virgin?
        freshnode = (lake_map[cnode] < 0) or (
            lake_map[cnode] >= grid.number_of_nodes)
        # Now, try to raise the level to this height, without running out
        # of water...
        z_increment = surf_z[nnode] - lake_water_level[cpit]
        (filled, z_increment) = _raise_water_level(cpit, cnode, z_increment,
                                                   flood_from_zero=freshnode, ...)
        lake_water_level[cpit] += z_increment
        if not filled:
            # Now, the lake can't grow any more! (...at the moment.)
            # there's still potential here for this lake to grow more, and
            # if it does, the first node to rise will be this one. So,
            # stick the node back in the local lake queue, along with all
            # the higher neighbors that never got raised to before we severed
            # the iteration:
            lake_q_dict[cpit].add_task(cnode, priority=lake_water_level[cpit])
            break

        # Now, a special case where we've filled the lake to this level,
        # but it turns out the next node is in fact the spill of an
        # adjacent lake...
        if 0 <= lake_map[nnode] < grid.number_of_nodes:
            # this condition only triggered (I think!) if we are now
            # raised to the spill level of an existing lake.
            # this then becomes break case (b) - except we cheat, and
            # don't actually break, since we're safe to continue raising
            # the new, composite lake as if it were a new one.
            cpit = _merge_two_lakes(
                this_pit=cpit, that_pit=lake_map[nnode])
            # we allow the loop to continue after this, provided the
            # queues and dicts are all correctly merged
            # in this case, leave the cnode as it is... they're all at the
            # same level now
            # allow a nghb check just in case we have a super funky geometry,
            # but very likely to find new neighbors
        cnode = nnode
        # ^Note, even in the merging case, we leave the actual sill as the next
        # node, so the current z_surf makes sense next time around

        # OK, so we've filled that node. Now let's consider where next...
        cnghbs = neighbors[cnode]
        # note that an "open" neighbor should potentially include the spill
        # node of an adjacent lake. This case dealt with in the raise step.
        for n in np.nditer(cnghbs):
            if not closednodes[n]:
                if lake_map[n] not in (cpit, pit_nghb_code, spill_code):
                    # ^virgin node (n==-1), or unclaimed nghb of other lake
                    # (0<=n<nnodes and n!=cpit), or future spill of another
                    # lake are all permitted.
                    # Now, if we've got to here, then all open neighbors must
                    # lead upwards; a spur rather than a sill will have already
                    # been explored from below. So if the topo falls, current
                    # node is a true sill
                    if surf_z[n] < surf_z[cnode]:
                        lake_map[cnode] = spill_code
                        # "magic" coding that lets this lake put "dibs" on that
                        # spill node.
##### is the spill coding strictly necessary?
                        lake_spill_node[cpit] = cnode
                        master_pit_q.add_task(
                            cpit, priority=lake_water_level[cpit])
                        tobreak = True
                        # ...but nevertheless, we want to keep loading the
                        # others so we can continue to use this list later if
                        # the lake acquires a new flux
                    else:
                        lake_q_dict[cpit].add_task(n, priority=surf_z[n])
                        lake_map[n] = pit_nghb_code


def _merge_two_lakes(grid, this_pit, that_pit, z_topo, lake_dict, lake_map,
                     lake_q_dict,
                     lake_water_level,
                     accum_area_at_pit,
                     vol_rem_at_pit,
                     lake_water_volume_balance,
                     accum_Ks_at_pit,
                     lake_is_full,
                     lake_spill_node):
    """
    Take two lakes that are known to be at the same level & in contact,
    and merge them.

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
    if z_topo[this_pit] < z_topo[that_pit]:
        low_pit = this_pit
        hi_pit = that_pit
    else:
        low_pit = that_pit
        hi_pit = this_pit

    # merge the queues
    lake_q_dict[low_pit].merge_queues(lake_q_dict[hi_pit])
    _ = lake_q_dict.pop(hi_pit)
    # merge the maps
    lake_map[lake_map == hi_pit] = low_pit
    # merge the various dicts
    for dct in (accum_area_at_pit, vol_rem_at_pit, lake_water_volume_balance,
                accum_Ks_at_pit):
        dct[low_pit] += dct.pop(hi_pit)
    for dct in (lake_water_level, lake_is_full, lake_spill_node):
        _ = dct.pop(hi_pit)
    lake_spill_node[low_pit] = -1  # just in case

    return low_pit


def _get_float_water_vol_balance_terms(
        water_vol_balance_terms, area_map, node):
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
    >>> type(c) is float
    True
    >>> type(K) is float
    True
    >>> c
    -3.
    >>> K
    0.

    >>> c, K = _get_float_water_vol_balance_terms((-1 * np.arange(4), 2.),
    ...                                           area_map, 3)
    >>> type(c) is float
    True
    >>> type(K) is float
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
                       accum_area_at_pit, vol_rem_at_pit, accum_Ks_at_pit,
                       lake_is_full, lake_spill_node,
                       flood_from_zero=True):
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
        The area of the cells at nodes across the grid.
    z_increment : float
        The total change in elevation between the current level and the target
        level.
    lake_map : array of ints
    water_vol_balance_terms : (c, K)

    accum_area_at_pit : dict
    vol_rem_at_pit : dict
    lake_is_full
    lake_spill_node

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
    >>> accum_Ks_at_pit = {5: 0.}
    >>> lake_is_full = {5: False}
    >>> lake_spill_node = {5: 8}
    >>> full, dz = _raise_water_level(cpit=5, cnode=3, area_map=area_map,
    ...                               z_increment=3., lake_map=lake_map,
    ...                               water_vol_balance_terms=(0., 0.),
    ...                               accum_area_at_pit=accum_area_at_pit,
    ...                               vol_rem_at_pit=vol_rem_at_pit,
    ...                               accum_Ks_at_pit=accum_Ks_at_pit,
    ...                               lake_is_full=lake_is_full,
    ...                               lake_spill_node=lake_spill_node,
    ...                               flood_from_zero=True)
    >>> full
    True
    >>> dz == 3.
    True
    >>> np.all_equal(lake_map, np.array([-1, -1, -1,
    ...                                   5,  5, -1,
    ...                                  -1, -1, -1]))
    True
    >>> np.isclose(accum_area_at_pit[5], 7.)
    True
    >>> np.isclose(vol_rem_at_pit[5], 79.)
    True
    >>> np.isclose(accum_Ks_at_pit[5], 0.)
    True
    >>> lake_is_full[5]
    False
    >>> lake_spill_node[5]
    8

    >>> full, dz = _raise_water_level(cpit=5, cnode=5, area_map=area_map,
    ...                               z_increment=4., lake_map=lake_map,
    ...                               water_vol_balance_terms=(-5., 0.),
    ...                               accum_area_at_pit=accum_area_at_pit,
    ...                               vol_rem_at_pit=vol_rem_at_pit,
    ...                               accum_Ks_at_pit=accum_Ks_at_pit,
    ...                               lake_is_full=lake_is_full,
    ...                               lake_spill_node=lake_spill_node,
    ...                               flood_from_zero=True)
    >>> full
    True
    >>> dz == 4.
    True
    >>> np.all_equal(lake_map, np.array([-1, -1, -1,
    ...                                   5,  5,  5,
    ...                                  -1, -1, -1]))
    True
    >>> np.isclose(accum_area_at_pit[5], 9.)
    True
    >>> np.isclose(vol_rem_at_pit[5], 33.)
    True
    >>> np.isclose(accum_Ks_at_pit[5], 0.)
    True
    >>> lake_is_full[5]
    False
    >>> lake_spill_node[5]
    8

    >>> full, dz = _raise_water_level(cpit=5, cnode=6, area_map=area_map,
    ...                               z_increment=1., lake_map=lake_map,
    ...                               water_vol_balance_terms=(-34., -3.),
    ...                               accum_area_at_pit=accum_area_at_pit,
    ...                               vol_rem_at_pit=vol_rem_at_pit,
    ...                               accum_Ks_at_pit=accum_Ks_at_pit,
    ...                               lake_is_full=lake_is_full,
    ...                               lake_spill_node=lake_spill_node,
    ...                               flood_from_zero=True)
    >>> full
    False
    >>> dz == 0.
    True
    >>> np.all_equal(lake_map, np.array([-1, -1, -1,
    ...                                   5,  5,  5,
    ...                                   5, -1, -1]))
    True
    >>> np.isclose(accum_area_at_pit[5], 10.)
    True
    >>> np.isclose(vol_rem_at_pit[5], -1.)
    True
    >>> np.isclose(accum_Ks_at_pit[5], 0.)
    True
    >>> lake_is_full[5]
    True
    >>> lake_spill_node[5]
    -1

    >>> lake_is_full[5] = False
    >>> lake_spill_node[5] = 8
    >>> vol_rem_at_pit[5] = 87.
    >>> full, dz = _raise_water_level(cpit=5, cnode=6, area_map=area_map,
    ...                               z_increment=2., lake_map=lake_map,
    ...                               water_vol_balance_terms=(-34., -3.),
    ...                               accum_area_at_pit=accum_area_at_pit,
    ...                               vol_rem_at_pit=vol_rem_at_pit,
    ...                               accum_Ks_at_pit=accum_Ks_at_pit,
    ...                               lake_is_full=lake_is_full,
    ...                               lake_spill_node=lake_spill_node,
    ...                               flood_from_zero=False)
    >>> full
    True
    >>> dz == 2.
    True
    >>> np.all_equal(lake_map, np.array([-1, -1, -1,
    ...                                   5,  5,  5,
    ...                                   5, -1, -1]))
    True
    >>> np.isclose(accum_area_at_pit[5], 10.)  # ...still
    True
    >>> np.isclose(vol_rem_at_pit[5], 61.)  # -6 loss, -20 rise
    True
    >>> np.isclose(accum_Ks_at_pit[5], -3.)
    True
    >>> lake_is_full[5]
    False
    >>> lake_spill_node[5]
    8

    >>> full, dz = _raise_water_level(cpit=5, cnode=2, area_map=area_map,
    ...                               z_increment=10., lake_map=lake_map,
    ...                               water_vol_balance_terms=(
    ...                                 np.arange(9), np.ones(9)),
    ...                               accum_area_at_pit=accum_area_at_pit,
    ...                               vol_rem_at_pit=vol_rem_at_pit,
    ...                               accum_Ks_at_pit=accum_Ks_at_pit,
    ...                               lake_is_full=lake_is_full,
    ...                               lake_spill_node=lake_spill_node,
    ...                               flood_from_zero=True)
    >>> full
    False
    >>> np.isclose(dz, 5.)

    >>> np.all_equal(lake_map, np.array([-1, -1,  5,
    ...                                   5,  5,  5,
    ...                                   5, -1, -1]))
    True
    >>> np.isclose(accum_area_at_pit[5], 12.)
    True
    >>> np.isclose(vol_rem_at_pit[5], 0.)  +4, then -5 (=-3+2) loss, -60 rise
    True
    >>> np.isclose(accum_Ks_at_pit[5], -3.)  # node not filled, so not changed
    True
    >>> lake_is_full[5]
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
            lake_is_full[cpit] = True
            lake_spill_node[cpit] = -1
            return (False, 0.)  # return a zero, but the node is still flooded

    V_gain_in_full_fill = (
        K_to_lift + accum_Ks_at_pit[cpit]) * z_increment  # likely <= 0.
    V_increment_to_fill -= V_gain_in_full_fill  # bigger if lossy
    # note that it's very possible that we fill the node then just add to
    # the total if the pit is, in fact, gaining... (careful not to double
    # account any gains wrt the flow routing)
    if vol_rem_at_pit[cpit] < V_increment_to_fill:
        # we don't have the water to get onto the next level
        frac = vol_rem_at_pit[cpit]/V_increment_to_fill
        lake_is_full[cpit] = True
        lake_spill_node[cpit] = -1
        z_available = frac * z_increment
        vol_rem_at_pit[cpit] = 0.
        return (False, z_available)
    else:
        # successful fill of node
        # increment the K term, so when we fill next this node is already
        # reflected:
        accum_Ks_at_pit[cpit] += K_to_lift
        vol_rem_at_pit[cpit] -= V_increment_to_fill
        return (True, z_increment)


def _route_outlet_to_next(outlet_ID, flow__receiver_nodes, z_surf, lake_map):
    """
    Take an outlet node, and then follow the steepest descent path until it
    reaches a distinct lake.
    """




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

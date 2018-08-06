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
                                        water_z, water_vol_balance_at_node,
                                        neighbors_at_nodes):
    """
    Take an outlet and its catchment nodes, then work outwards from the
    pit nodes using the volumes of water at each to work out which nodes
    are connected to which, and where is inundated.

    Assumes we already know that the pit doesn't overtop, i.e., is its own
    little endorheic basin.

    DO NOT INCLINE THE WATER SURFACES!!

    Idea here is that water_vol_balance_at_node IS A FUNC OR ARRAY, such that
    it can be sensitive to depth if necessary.
    Provide water_vol_balance_at_node as [c, K], where balance = K*depth + c.
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
    accum_ks_at_pit = {}
    accum_cs_at_pit = {}
    lake_is_full = {}
    for pit, vol, A, level in zip(pit_nodes_in_order, Vw_at_pits_in_order,
                                  area_map[pit_node_sort],
                                  zwater_nodes_in_order):
        vol_rem_at_pit[pit] = vol
        accum_area_at_pit[pit] = A
        lake_water_level[pit] = level
        lake_is_full[pit] = False
        lake_water_volume_balance[pit] = 0.
        accum_ks_at_pit[pit] = 0.
        accum_cs_at_pit[pit] = 0.
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
                accum_ks_at_pit
                accum_cs_at_pit
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
                    accum_ks_at_pit
                    accum_cs_at_pit
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
                                                 water_vol_balance_at_node,
                                                 accum_ks_at_pit, c_to_add)
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
                        accum_ks_at_pit[currentnode_pit] * z_increment)
            # now, simply move to this next lowest node
            current_node = nextlowest


def _raise_lake_to_limit(grid, current_pit, lake_dict, lake_map, master_pit_q,
                         lake_q_dict, init_water_surface,
                         lake_water_level,
                         accum_area_at_pit,
                         vol_rem_at_pit,
                         lake_water_volume_balance,
                         accum_ks_at_pit,
                         accum_cs_at_pit,
                         lake_is_full,
                         lake_spill_node,
                         water_vol_balance_at_node,
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

    Parameters
    ----------
    neighbors : (nnodes, max_nneighbors) array
        The neighbors at each node on the grid.
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
        A dict of queues holding nodes known to be on the perimeter of each lake.
        We hold this information since we return to each lake repeatedly to
        update the elevs. If the lake hasn't been filled before, this will
        just contain the pit node.
    init_water_surface : array
        The original topo of the water surface at the start of the step,
        i.e., before we start moving the lake water levels this time around.
    lake_water_level : dict
        The current (transient) water level at each lake.
     accum_area_at_pit,
     vol_rem_at_pit,
     lake_water_volume_balance,
     accum_ks_at_pit,
     accum_cs_at_pit,
     lake_is_full,

     lake_spill_node : dict
        the current spill node of each lake. -1 if not defined.
    water_vol_balance_at_node : (float, float)
        polyparams for a linear fit of the water vol balance func.
    """
    # We work upwards in elev from the current level, raising the level to
    # the next lowest node. We are looking for the first sign of flow
    # going into a node that isn't already in the queue that is *downhill*
    # of the current elevation. The idea here is that either a node is
    # accessible from below via another route, in which case it's already in
    # the queue, or it's over a saddle, and thus belongs to a separate pit.
    tobreak = False  # delayed break flag
    cpit = current_pit  # shorthand
#     abs_elev_limit = master_pit_q.peek_at_task()
# #### probably redundant
    pit_nghb_code = cpit + grid.number_of_nodes
    # ^this code used to flag possible nodes to inundate, but that aren't
    # flooded yet.
    surf_z = init_water_surface
    while not tobreak:
        cnode = lake_q_dict[cpit].pop_task()
        # note that this is the next node we will try to flood, not the one
        # we are currently flooding. But conversely, a node that has the
        # water at its level does count as flooded by this lake (i.e., a
        # lake sill will be incorporated into that lake once the water
        # reaches it)
        # Now, try to raise the level to this height, without running out
        # of water...
        z_increment = surf_z[cnode] - lake_water_level[cpit]
        # # Check we aren't about to attempt a raise over and above the
        # # absolute limit:
        # if surf_z[cnode] > abs_elev_limit:
        #     z_increment = abs_elev_limit - lake_water_level[cpit]
        #     master_pit_q.add_task(cpit, priority=lake_water_level[cpit])
        #     # put level back in the list
        #     tobreak = True
        V_increment_to_fill = (
            z_increment * accum_area_at_pit[cpit] -
            lake_water_volume_balance[cpit])
        # & remember to update the vol_balance below
        if V_increment_to_fill < 0.:
            # this really should not happen!!
            raise AssertionError("lake is gaining too much direct input!!")
        # Note that the addition of the vol balance makes it harder to
        # actually fill the step in most cases.

        if vol_rem_at_pit[cpit] <= V_increment_to_fill:
            # we don't have the water to get onto the next level
            lake_is_full[cpit] = True
            lake_spill_node[cpit] = -1
            z_available = vol_rem_at_pit[cpit]/accum_area_at_pit[cpit]
            lake_water_level[cpit] += z_available
            # no need to update water balance terms, since (a) this lake is
            # done, and (b) k's anc c's don't change anyway  as lake area is
            # static
            # Now, the lake can't grow any more! (...at the moment.)
            # there's still potential here for this lake to grow more, and
            # if it does, the first node to rise will be this one. So,
            # stick the node back in the local lake queue, along with all
            # the higher neighbors that never got raised to before we severed
            # the iteration:
            lake_q_dict[cpit].add_task(cnode, priority=lake_water_level[cpit])
            # But do not re-list this in the main queue!
            # Update the loss due to the increased depth:
            lake_water_volume_balance[cpit] += (
                accum_ks_at_pit[c_pit] * z_available)
            # Then simply...
            break  # break cond (a)
        else:
            # successful fill of node
            vol_rem_at_pit[cpit] -= V_increment_to_fill
            lake_water_level[cpit] += z_increment
            accum_area_at_pit[cpit] += grid.cell_area_at_node[cnode]

            # Now, a special case where we've filled the lake to this level,
            # but it turns out this node is in fact the spill of an
            # adjacent lake...
            if 0 <= lake_map[cnode] < grid.number_of_nodes:
                # this condition only triggered (I think!) if we are now
                # raised to the spill level of an existing lake.
                # this then becomes break case (b) - except we cheat, and
                # don't actually break, since we're safe to continue raising
                # the new, composite lake as if it were a new one.
                _some_func_to_merge_lakes(this_pit, that_pit)
                # we allow the loop to continue after this, provided the
                # queues and dicts are all correctly merged
            else:  # we only want to add the node if it's fresh!
                # update the dicts describing the lake, as it floods onto
                # the new surface.
                _propagate_lake_to_next_node(cpit, cnode,
                                             grid.cell_area_at_node,
                                             z_increment,
                                             accum_area_at_pit,
                                             lake_water_volume_balance,
                                             water_vol_balance_at_node,
                                             accum_ks_at_pit, c_to_add)
                # incorporate the node into the lake:
                lake_map[cnode] = cpit
                # (this redundant if we merged lakes)

        # OK, so we've filled that node. Now let's consider where next...
        cnghbs = neighbors[cnode]
        # note that an "open" neighbor should potentially include the spill
        # node of an adjacent lake. This case dealt with in the raise step.
        for n in np.nditer(cnghbs):
            if lake_map[n] not in (cpit, pit_nghb_code):
                # ^virgin node (n==-1), or spill of other lake (0<=n<nnodes and
                # n!=cpit), or node is already in nghb queue of another lake,
                # but not actually flooded (n>=nnodes and n!=pit_nghb_code).
                # Now, if we've got to here, then all open neighbors must lead
                # upwards; a spur rather than a sill will have already been
                # explored from below. So if the topo falls, that's a true sill
                if surf_z[n] < surf_z[cnode]:
                    lake_spill_node[cpit] = cnode
                    master_pit_q.add_task(
                        cpit, priority=lake_water_level[cpit])
                    tobreak = True
                    # ...but nevertheless, we want to keep loading the others
                    # so we can continue to use this list later if the lake
                    # acquires a new flux
                else:
                    lake_q_dict[cpit].add_task(n, priority=surf_z[n])
                    lake_map[n] = pit_nghb_code
                    # we're perfectly happy to add nghbs that exceed the
                    # (current) abs_elev_limit, since we are likely to be
                    # back in this lake at some point


def _merge_two_lakes(grid, this_pit, that_pit, lake_dict, lake_map,
                     lake_q_dict,
                     lake_water_level,
                     accum_area_at_pit,
                     vol_rem_at_pit,
                     lake_water_volume_balance,
                     accum_ks_at_pit,
                     accum_cs_at_pit,
                     lake_is_full,
                     lake_spill_node):
    """
    Take two lakes that are known to be at the same level & in contact,
    and merge them.
    """
    # Note we *will* still have stuff in the master_pit_q that has the ID of
    # a replaced lake. Deal with this by spotting that those lake codes no
    # longer work as dict keys...
    # Check they've got to the same level (really ought to have!)
    assert np.isclose(lake_water_level[this_pit], lake_water_level[that_pit])
    _ = lake_water_level.pop(that_pit)
    # merge the queues
    lake_q_dict[this_pit].merge_queues(lake_q_dict[that_pit])
    _ = lake_q_dict.pop(that_pit)
    # merge the maps
    lake_map[lake_map == that_pit] = this_pit
    # merge the various dicts
    for dct in (accum_area_at_pit, vol_rem_at_pit, lake_water_volume_balance,
                accum_ks_at_pit):
        dct[this_pit] += dct.pop(that_pit)
    # the c dict is OK (not adding a new node); just remove it in that_pit
    for dct in (accum_cs_at_pit, lake_is_full, lake_spill_node):
        _ = dct.pop(that_pit)
    lake_spill_node[this_pit] = -1


def _propagate_lake_to_next_node(current_pit_id, node_to_flood,
                                 area_map, z_increment,
                                 accum_area_at_pit, lake_water_volume_balance,
                                 water_vol_balance_at_node,
                                 accum_ks_at_pit, c_to_add):
    """
    Helper fn to allow water to flood onto an adjacent node within a filling
    lake.
    """
    accum_area_at_pit[current_pit_id] += area_map[node_to_flood]
    # Update the water balance terms as we spill over
    if type(water_vol_balance_at_node[1]) is np.array:
        accum_ks_at_pit[current_pit_id] += (
            water_vol_balance_at_node[1][node_to_flood])
    else:
        accum_ks_at_pit[current_pit_id] += (
            water_vol_balance_at_node[1])
    if type(water_vol_balance_at_node[0]) is np.array:
        c_to_add = water_vol_balance_at_node[0][node_to_flood]
    else:
        c_to_add = water_vol_balance_at_node[0]
    lake_water_volume_balance[current_pit_id] += (
        accum_ks_at_pit[current_pit_id] * z_increment + c_to_add)

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

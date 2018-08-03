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

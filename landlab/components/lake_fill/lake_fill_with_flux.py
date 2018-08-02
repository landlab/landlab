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
    """
    disch_map = mg.zeros('node', dtype=float)  # really a vol map!
    lake_map = mg.zeros('node', dtype=int)
    lake_map.fill(-1)  # this will store the unique lake IDs
    area_map = mg.cell_area_at_node  # we will update this as we go...
    # sort the pit nodes. We need to start at the lowest.
    zsurf_pit_nodes = surface_z[pit_nodes]
    pit_node_sort = np.argsort(zsurf_pit_nodes)
    pit_nodes_in_order = pit_nodes[pit_node_sort]
    Vw_at_pits_in_order = Vw_at_pits[pit_node_sort]
    vol_rem_at_pit = {}
    accum_area_at_pit = {}
    lake_water_level = {}
    lake_is_full = {}
    for pit, vol, A, level in zip(pit_nodes_in_order, Vw_at_pits_in_order,
                                  area_map[pit_node_sort],
                                  zsurf_pit_nodes[pit_node_sort]):
        vol_rem_at_pit[pit] = vol
        accum_area_at_pit[pit] = A
        lake_water_level[pit] = level
        lake_is_full[pit] = False
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
        flood_order.add_task(pit, priority=surface_z[pit])

    current_node = flood_order.pop_task()
    while True:
        currentnode_pit = lake_map[current_node]
        if lake_is_full[currentnode_pit]:
            # stop treating this pit entirely; just whip through any nodes
            # on its perimeter w/o action
            try:
                current_node = flood_order.pop_task()
            except KeyError:
                break
            continue
        # for the purposes of the logic, a flooded node (i.e., tagged w ID)
        # is one which has AT LEAST SOME inundation. Does NOT have to be
        # fully filled
        zsurf_node = surface_z[current_node]
        # at this point we know the water is lapping at the foot of this node
        node_neighbors = neighbors_at_nodes[current_node]
        currentnode_A = accum_area_at_pit[currentnode_pit]  # already incremented; do this!
        currentnode_V_rem = vol_rem_at_pit[currentnode_pit]  # already incremented

        # load up the queue. Likely one of these will be next, so do it 1st
        for nghb in np.nditer(node_neighbors):
            if lake_map[nghb] < 0:
                zsurf_nghb = surface_z[nghb]
                if zsurf_nghb >= zsurf_node:
                    flood_order.add_task(nghb, priority=zsurf_nghb)
                    lake_map[nghb] = currentnode_pit
            else:
                pass ######FOUND A SILL, TREAT HERE...? -> YES

        # now, pull the lowest thing out of the queue:
        try:
            nextlowest = flood_order.pop_task()  # TREAT THE ERROR
        except KeyError:
            break  # filled the pit completely
        if lake_map[nghb] >= 0:  # no risk of finding itself; it's above level
            # found an existing flooded surface, i.e., A SILL
            # Note possible for this to be a sill now when it wasn't before...
            # COME BACK AND SAVE THIS INFO
            # flag the sill, then route the discharge down into the next
            # depression?? (grab the flowdirs?)
            # do we mean continue? maybe still need to flood up. Or do later?
            # THINK THIS IS REDUNDANT SO
            raise AssertionError
        zsurf_next = surface_z[nextlowest]
        z_increment = zsurf_next - zsurf_node
        # try to flood up to next
        V_increment_to_fill = z_increment * currentnode_A

        if currentnode_V_rem <= V_increment_to_fill:
            # we don't have the water to get onto the next level
            lake_is_full[currentnode_pit] = True
            z_available = currentnode_V_rem/currentnode_A
            lake_water_level[currentnode_pit] += z_available
        else:
            # successful fill of node
            vol_rem_at_pit[currentnode_pit] -= V_increment_to_fill
            lake_water_level[currentnode_pit] += z_increment
            # propagate the lake info into this next node, if it's
            # contiguous...
            if lake_map[current_node] == lake_map[nextlowest]:
                accum_area_at_pit[lake_ID] += area_map[nextlowest]
            else:
                # NOT contiguous. In this case, we put it back in the
                # queue, as it has more growing to do to reach an actual
                # neighbor:
                flood_order.add_task(current_node, priority=zsurf_node)
                # note that due to the ordered queue, we will handle the
                # two pits that must be dealt with here in an alternating
                # fashion. This is important, since we are bringing both
                # up to the same level simulataneously.
        # now, simply move to this next lowest node
        current_node = nextlowest




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

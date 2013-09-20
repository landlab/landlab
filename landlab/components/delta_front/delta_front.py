import numpy as np
import numpy.ma as ma

"""
Preliminary work for a submarine delta foreset sed router. Not yet functioning!
DEJH, Sept 2013
"""

class Littoral(object):
    def __init__(self, grid, dry_land_Qs_array, elev):
        #dry_land_Qs_array is a np array of -1s for any cell <SL, with Qs for any dry land cell.
        #Very important that when sed is routed on dry land, sea cells act as hard boundaries, BC=4.
        self.dry_land_Qs = dry_land_Qs_array[dry_land_Qs_array>0]
        
        
    def assign_Qs_to_littoral_cells(self, grid, dry_land_Qs_array, elev):
        self.dry_land_Qs = dry_land_Qs_array[dry_land_Qs_array>0]
        cell_neighbors = grid.get_neighbor_list(self.dry_land_Qs)
        cell_neighbor_heights = elev[cell_neighbors]
        number_wet_cell_neighbors = np.sum(cell_neighbor_heights<0, axis=1)
        wet_neighbor_nonzero = np.where(number_wet_cell_neighbors>0)
        wet_neighbor_count_nonzero = number_wet_cell_neighbors[wet_neighbor_nonzero]
        Qs_to_each_neighbor = dry_land_Qs_array[wet_neighbor_nonzero]/wet_neighbor_count_nonzero
        wet_cell_neighbors = np.where(cell_neighbor_heights<0,cell_neighbors,-1)[wet_neighbor_nonzero]
        self.Qs_in_source_cell_extended = np.zeros(len(dry_land_Qs_array)+1)
        self.Qs_in_source_cell_extended[wet_cell_neighbors] += Qs_to_each_neighbor #broadcasting should sort this out automatically
        self.Qs_in_source_cell = self.Qs_in_source_cell_extended[:-1]
        
        return self.Qs_in_source_cell
        #This returned array is filled with zeros, except where sed has been discharged directly into the sea cells, where it
        #is equal to the total sed flux entering that cell.

class ForesetAggrade(object):
    def __init__(self, grid, source_cells_Qs, elev):
        self.tan_repose_angle = np.tan(32.*np.pi/180)
        #source_cells_Qs is a np array of 0s except where a (submarine) source cell, when it contains the Qs delivered to that cell from land.
        
    def aggrade_front(self, grid, tstep, source_cells_Qs, elev, SL): #ensure Qs and tstep units match!
        
        self.total_sed_supplied_in_tstep = source_cells_Qs*tstep
        self.Qs_sort_order = np.argsort(source_cells_Qs)[::-1] #descending order
        self.Qs_sort_order = self.Qs_sort_order[:np.count_nonzero(self.Qs_sort_order>0)]
        for i in self.Qs_sort_order:
            subaerial_nodes = elev>=SL
            subsurface_elev_array = ma.array(elev, subaerial_nodes)
            xy_tuple = (grid.node_x[i], grid.node_y[i])
            distance_map = grid.get_distances_of_nodes_to_point(xy_tuple)
            loop_number = 0
            closest_node_list = ma.argsort(ma.masked_array(distance_map, mask=subsurface_elev_array.mask))
            smooth_cone_elev_from_apex = subsurface_elev_array[i]-distance_map*self.tan_repose_angle
            while 1:
                filled_all_cells_flag = 0
                accom_space_at_controlling_node = SL - subsurface_elev_array[closest_node_list[loop_number]]
                new_max_cone_surface_elev = smooth_cone_elev_from_apex + accom_space_at_controlling_node
                subsurface_elev_array.mask = (elev>=SL or new_max_cone_surface_elev<elev)
                depth_of_accom_space = new_max_cone_surface_elev - subsurface_elev_array
                accom_depth_order = ma.argsort(depth_of_accom_space)[::-1]
                #Vectorised method to calc fill volumes:
                area_to_fill = ma.cumsum(grid.cellarea[accom_depth_order])
                differential_depths = ma.empty_like(depth_of_accom_space)
                differential_depths[:-1] = depth_of_accom_space[accom_depth_order[:-1]] - depth_of_accom_space[accom_depth_order[1:]]
                differential_depths[-1] = depth_of_accom_space[accom_depth_order[-1]]
                incremental_volumes = ma.cumsum(differential_depths*area_to_fill)
                match_position_of_Qs_in = ma.searchsorted(incremental_volumes, self.total_sed_supplied_in_tstep[i])
                try:
                    depths_to_add = depth_of_accom_space-depth_of_accom_space[match_position_of_Qs_in]
                except:
                    depths_to_add = depth_of_accom_space-depth_of_accom_space[match_position_of_Qs_in-1]
                    filled_all_cells_flag = 1
                depths_to_add = depths_to_add[ma.where(depths_to_add>=0)]
                if not filled_all_cells_flag:
                    depths_to_add += (self.total_sed_supplied_in_tstep[i] - incremental_volumes[match_position_of_Qs_in-1])/area_to_fill[match_position_of_Qs_in-1]
                    subsurface_elev_array[accom_depth_order[len(depths_to_add)]] = depths_to_add
                    self.total_sed_supplied_in_tstep[i] = 0
                    break
                else:
                    subsurface_elev_array[accom_depth_order] = depths_to_add
                    self.total_sed_supplied_in_tstep[i] -= incremental_volumes[-1]
                    loop_number += 1
        
        return elev

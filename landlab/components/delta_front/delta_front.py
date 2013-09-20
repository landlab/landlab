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
        #source_cells_Qs is a np array of -1s except where a source cell, when it contains the Qs in
        
    def aggrade_front(self, grid, tstep, source_cells_Qs, elev, SL): #ensure Qs and tstep units match!
        
        self.total_sed_supplied_in_tstep = source_cells_Qs*tstep
        self.Qs_sort_order = np.argsort(source_cells_Qs)[::-1] #descending order
        self.Qs_sort_order = self.Qs_sort_order[:np.count_nonzero(self.Qs_sort_order>0)]
        for i in self.Qs_sort_order:
            subaerial_nodes = elev>SL
            subsurface_elev_array = ma.array(elev, subaerial_nodes)
            xy_tuple = (grid.node_x[i], grid.node_y[i])
            distance_map = grid.get_distances_of_nodes_to_point(xy_tuple)
            #cone_elev_from_current_source = subsurface_elev_array[i]-grid.get_distances_of_nodes_to_point(xy_tuple)*self.tan_repose_angle
            max_cone_elev_to_surface = -distance_map*self.tan_repose_angle
            subsurface_elev_array[max_cone_elev_to_surface<subsurface_elev_array] = ma.masked
            #the subsurface_elev_array is now masked to the cells into which sed can be added to form the cone
            depth_of_accom_space = max_cone_elev_to_surface - subsurface_elev_array
            accom_depth_order = ma.argsort(depth_of_accom_space)[::-1]
            for j in range(ma.count_nonzero(depth_of_accom_space)):
                #this try loop will break when the final 
                try:
                    height_to_raise = depth_of_accom_space[accom_depth_order[j]]-depth_of_accom_space[accom_depth_order[j+1]]
                except:
                    break
                else:
                    area_raised = np.sum(grid.cellarea[accom_depth_order[:j]])
                    if height_to_raise*area_raised < self.total_sed_supplied_in_tstep[i]:
                        self.total_sed_supplied_in_tstep[i] -= height_to_raise*area_raised
                        subsurface_elev_array[accom_depth_order[:j]] += height_to_raise
                    else:
                        height_to_raise = self.total_sed_supplied_in_tstep[i]/area_raised
                        subsurface_elev_array[accom_depth_order[:j]] += height_to_raise
                        break
            if self.total_sed_supplied_in_tstep[i] > 0:
                if np.sum(grid.cellarea[accom_depth_order])*depth_of_accom_space[accom_depth_order[-1]] > self.total_sed_supplied_in_tstep[i]:
                    subsurface_elev_array += self.total_sed_supplied_in_tstep[i]/np.sum(grid.cellarea[accom_depth_order])
                else:
                    subsurface_elev_array += depth_of_accom_space[accom_depth_order[-1]]
                    self.total_sed_supplied_in_tstep[i] -= np.sum(grid.cellarea[accom_depth_order])*depth_of_accom_space[accom_depth_order[-1]]
                    #Then need to do further raising of next cells out... but at least it's smooth in the masked grid already.
                    loop_number = 1 #1 is right; don't want to select the apex node itself!
                    closest_node_list = ma.argsort(ma.masked_array(distance_map, mask=subsurface_elev_array.mask))
                    while 1:
                        accom_space_controlling_node = SL - subsurface_elev_array[closest_node_list[loop_number]]
                        new_max_cone_surface_elev = max_cone_elev_to_surface + accom_space_controlling_node
                        subsurface_elev_array.mask = (elev>=SL or new_max_cone_surface_elev<elev)
                        #Now more or less same process as above:
                        depth_of_accom_space = new_max_cone_surface_elev - subsurface_elev_array
                        accom_depth_order = ma.argsort(depth_of_accom_space)[::-1]
                        for k in range(ma.count_nonzero(depth_of_accom_space)):
                            try:
                                height_to_raise = depth_of_accom_space[accom_depth_order[k]]-depth_of_accom_space[accom_depth_order[k+1]]
                            except:
                                break
                            else:
                                area_raised = np.sum(grid.cellarea[accom_depth_order[:k]])
                                if height_to_raise*area_raised < self.total_sed_supplied_in_tstep[i]:
                                    self.total_sed_supplied_in_tstep[i] -= height_to_raise*area_raised
                                    subsurface_elev_array[accom_depth_order[:k]] += height_to_raise
                                else:
                                    height_to_raise = self.total_sed_supplied_in_tstep[i]/area_raised
                                    subsurface_elev_array[accom_depth_order[:k]] += height_to_raise
                                    break
                        if self.total_sed_supplied_in_tstep[i] > 0:
                            if np.sum(grid.cellarea[accom_depth_order])*depth_of_accom_space[accom_depth_order[-1]] > self.total_sed_supplied_in_tstep[i]:
                                subsurface_elev_array += self.total_sed_supplied_in_tstep[i]/np.sum(grid.cellarea[accom_depth_order])
                                break
                            else:
                                subsurface_elev_array += depth_of_accom_space[accom_depth_order[-1]]
                                self.total_sed_supplied_in_tstep[i] -= np.sum(grid.cellarea[accom_depth_order])*depth_of_accom_space[accom_depth_order[-1]]
                                loop_number += 1
                        else:
                            break
                        #We can almost certainly combine the first loop into this 2nd while loop
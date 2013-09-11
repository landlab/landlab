import numpy as np

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
        Qs_in_source_cell_extended[wet_cell_neighbors] += Qs_to_each_neighbor #broadcasting should sort this out automatically
        self.Qs_in_source_cell = Qs_in_source_cell_extended[:-1]
        
        return self.Qs_in_source_cell
        #This returned array is filled with zeros, except where sed has been discharged directly into the sea cells, where it
        #is equal to the total sed flux entering that cell.

class ForesetAggrade(object):
    def __init__(self, grid, source_cells_Qs, elev):
        self.repose_angle = np.tan(32.*np.pi/180)
        #source_cells_Qs is a np array of -1s except where a source cell, when it contains the Qs in
        
    def aggrade_front(self, grid, source_cells_Qs, elev):
        
        self.Qs_sort_order = np.argsort(source_cells_Qs)
        self.Qs_sort_order = self.Qs_sort_order[:np.count_nonzero(Qs_sort_order>0)]
        
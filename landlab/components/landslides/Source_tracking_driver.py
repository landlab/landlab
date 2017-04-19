# -*- coding: utf-8 -*-

"""
Created on Mon Apr 17 22:47:22 2017

@author: saisiddu
"""
from landlab.io.esri_ascii import read_esri_ascii
import SourceTrackingAlgorithm as STA
import cPickle as pickle


# %%
# Import input files. Create a RasterModelGrid object 'grid' and assign
# input data as fields or assigned to a variable
grid, z = read_esri_ascii('./Input_files/elevation.txt',
                          name='topographic__elevation')
grid.set_nodata_nodes_to_closed(grid['node']['topographic__elevation'], -9999.)
grid, flow_dir_arc = read_esri_ascii('./Input_files/flow_direction.txt',
                                     name='flow_dir', grid=grid)
grid, hsd_ids = read_esri_ascii('./Input_files/vic_idsnoca.txt',
                                name='hsd_id', grid=grid)
grid, slp_g16 = read_esri_ascii('./Input_files/slp_g16msk.txt',
                                name='slp_g16', grid=grid)
hsd_ids = hsd_ids.astype(int)
# Define model domain
grid.set_closed_boundaries_at_grid_edges(True, True, True, True)
grid.set_nodata_nodes_to_closed(grid['node']['flow_dir'], -9999.)
grid.set_nodata_nodes_to_closed(grid['node']['slp_g16'], -9999.)


# %%
r = STA.convert_arc_flow_directions_to_landlab_node_ids(grid, flow_dir_arc)
(hsd_upstr, flow_accum) = STA.track_source(grid, r, hsd_ids)
(uniq_ids, coeff) = STA.find_unique_upstream_hsd_ids_and_fractions(hsd_upstr)


# %%
# Plot or/and Save files
# Saving dict{MD id: unique HSD ids}
pickle.dump(uniq_ids, open("dict_uniq_ids.p", "wb"))
# Saving dict{MD id: Fraction of unique HSD ids
# contributing to this MD id}
pickle.dump(coeff, open("dict_coeff.p", "wb"))

# -*- coding: utf-8 -*-
"""
This code outlines very a very basic use case for the NetworkSedimentTransporter
component. 

Created on Sun May 20 15:54:03 2018

@authors: Jon Czuba, Allison Pfeiffer, Katy Barnhart
"""
import numpy as np
from landlab.grid.network import NetworkModelGrid
from landlab.item_collection import ItemCollection
from landlab.plot import graph

# %% Set the geometry using Network model grid (should be able to read in a shapefile here)

y_of_node = (0, 1, 2, 2, 3, 4, 4, 1.25)
x_of_node = (0, 0, 1, -0.5, -1, 0.5, -1.5, -1)

nodes_at_link = ((1, 0), (2, 1), (1, 7), (3, 1), (3, 4), (4, 5), (4, 6))

grid = NetworkModelGrid((y_of_node, x_of_node), nodes_at_link)
graph.plot_graph(grid, at='node,link')

grid.at_node['topographic__elevation'] = [0., 1., 3., 2., 3., 4., 4.1, 5.]
#^ in order for the FlowDirector and FlowAccumulator to work properly with the network model grid
# I had to change the last two elements from 4. and 2. to 4.1 and 5. in order to avoid slopes <=0

area = grid.add_ones('cell_area_at_node', at = 'node')


#%% Set geometry for each link

# Ultimately, map between flow accumulator and shapefile reader info...
# map_upstream_node_to_link

grid.at_link['drainage_area'] = [100e+6,10e+6,70e+6,20e+6,70e+6,30e+6,40e+6] # m
grid.at_link['channel_slope'] = [0.01,0.02,0.01,0.02,0.02,0.03,0.03]  
grid.at_link['link_length'] = [100,100,100,100,100,100,100] # m

grid.at_link['channel_width'] = 15 * np.ones(np.size(grid.at_link['drainage_area'])) # m REPLACE with something hydraulically meaningful
grid.at_link['channel_width'][3]= 10
# modify elevations so they are consistent with adjusted slopes

## Basic parameters

g = 9.81 # m/s2
rho = 1000 # kg/m3
theta = 0.5

Lp = 0.3  #porosity of the bed material

# %% initialize bed sediment (will become its own component)

# Ultimately,
# parcels = SedimentParcels(grid,initialization_info_including_future_forcing)
    # parcels have 'time added' and 'starting location/link' as an attribute. 

element_id = np.array([0, 0, 1, 1, 2, 2, 3, 3, 4, 4, 5, 5, 6, 6, 1], dtype=int) # current link for each parcel
starting_link = element_id # starting link for each parcel
#parcel_id = np.arange(0,np.size(element_id)) # each parcel has a unique identifier
np.random.seed(0)
time_arrival_in_link = np.random.rand(np.size(element_id)) # time of arrival in each link -- larger numbers are younger
volume = np.ones(np.size(element_id)) # (m3) the volume of each parcel
D = 0.05 * np.ones(np.size(element_id)) # (m) the diameter of grains in each parcel
lithology = ['quartzite']*np.size(element_id) # a lithology descriptor for each parcel
abrasion_rate = 0 * np.ones(np.size(element_id)) # 0 = no abrasion; abrasion rates are positive mass loss coefficients
active_layer = np.ones(np.size(element_id)) # 1 = active/surface layer; 0 = subsurface layer
density = 2650 * np.ones(np.size(element_id)) # (kg/m3) 

location_in_link = np.zeros(np.size(element_id)) # [0 1], 0 is upstream end of link, 1 is downstream end

D[0] = 0.075
D[5] = 0.0001 # make one of them sand

volume[2] = 0.3

data = {'starting_link': starting_link,
        'volume': volume,
        'D': D,
        'lithology': lithology,
        'time_arrival_in_link': time_arrival_in_link,
        'active_layer': active_layer,
        'density': density,
        'location_in_link': location_in_link
        'abrasion_rate': abrasion_rate}
        
parcels = ItemCollection(grid,
    data = data,
    grid_element ='link',
    element_id = element_id)

# Add parcels in at a given time --> attribute in the item collection


# %% Flow parameters (this will happen in the sq.run_one_step and sc.run_one_step)

timesteps = 1;

Qgage = 200. # for Methow, this is ~4*mean flow
dt = 60*60*24; # (seconds) daily timestep
Bgage=30.906*Qgage**0.1215; # (m)
Hgage=0.2703*Qgage**0.3447; # (m)
Agage=4.5895e+9;# (m2)

B=((np.tile(Bgage,(grid.number_of_links))/(Agage**0.5))
  *np.tile(grid.at_link['drainage_area'],(timesteps))**0.5)

H=((np.tile(Hgage,(grid.number_of_links))/(Agage**0.4))
  *np.tile(grid.at_link['drainage_area'],(timesteps))**0.4)


Btmax=np.amax(B, axis = 0)


# %% Instantiate component(s)
#dis = ExteralDischargeSetter(grid, filename, model='dhsvm')
        
#sq = SyntheticDischargeMaker(discharge,drainage_area) # OR read DHSVM. 
    # define a surface_water_discharge for each link and each timestep
                    
#sc = SyntheticChannelGeomMaker(hydraulic_geometry_scaling_rules,discharge)
    # 

nst = NetworkSedimentTransporter(grid,parcels,depth...?)

# %% Run the component(s)

timesteps = 10;

for t in range(timesteps):
    
    # move any sediment additions from forcing Item collector to bed item collector
    
    # sq.run_one_step  
    #   will assign discharge for each reach

    
    # sc.run_one_step   
    #   will assign flow depth for each reach (for this timestep)

    
    # Run our component
   nst.run_one_step

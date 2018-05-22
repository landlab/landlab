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

grid.at_node['topographic__elevation'] = [0., 1., 3., 2., 3., 4., 4., 2.]

area = grid.add_ones('cell_area_at_node', at = 'node')


#%% Set geometry for each link

# Ultimately, map between flow accumulator and shapefile reader info...
# map_upstream_node_to_link

grid.at_link['drainage_area'] = [1000,100,700,200,700,300,400] # m
grid.at_link['channel_slope'] = [0.01,0.02,0.01,0.02,0.02,0.03,0.03]  
grid.at_link['link_length'] = [100,100,100,100,100,100,100] # m

# modify elevations so they are consistent with adjusted slopes


# %% initialize bed sediment (will become its own component)

# Ultimately,
# parcels = SedimentParcels(grid,initialization_info_including_future_forcing)
    # parcels have 'time added' and 'starting location/link' as an attribute. 

element_id = np.array([0, 0, 1, 1, 2, 2, 3, 3, 4, 4, 5, 5, 6, 6], dtype=int) # current link for each parcel
starting_link = element_id # starting link for each parcel
time_arrival_in_link = np.zeros(np.size(element_id)) # time of arrival in each link
volume = np.ones(np.size(element_id)) # (m3) the volume of each parcel
D = 0.05 * np.ones(np.size(element_id)) # (m) the diameter of grains in each parcel
lithology = ['quartzite']*np.size(element_id) # a lithology descriptor for each parcel
active_layer = np.ones(np.size(element_id)) # 1 = active/surface layer; 0 = subsurface layer

data = {'starting_link': starting_link,
        'volume': volume,
        'D': D,
        'lithology': lithology,
        'time_arrival_in_link': time_arrival_in_link,
        'active_layer': active_layer}
        
parcels = ItemCollection(grid,
    data = data,
    grid_element ='link',
    element_id = element_id)


# Add parcels in at a given time --> attribute in the item collection

# %% Instantiate component(s)
#dis = ExteralDischargeSetter(grid, filename, model='dhsvm')
        
#sq = SyntheticDischargeMaker(discharge,drainage_area) # OR read DHSVM. 
    # define a surface_water_discharge for each link and each timestep
                    
#sc = SyntheticChannelGeomMaker(hydraulic_geometry_scaling_rules,discharge)
    # 

nst = NetworkSedimentTransporter(so many inputs!)

# %% Run the component(s)

for t in range(timesteps):
    
    # move any sediment additions from forcing Item collector to bed item collector
    
    # sq.run_one_step  
    Qgage = 200. # for Methow, this is ~4*mean flow
    dt = 60*60*24; # (seconds) daily timestep
    
    # sc.run_one_step   
    Bgage=30.906.*Qgage.^0.1215;
    Hgage=0.2703.*Qgage.^0.3447;
    Agage=4.5895e+9;#m2
    B=(repmat(Bgage,1,LinkNum)./(Agage.^0.5)).*repmat(usarea',timesteps,1).^0.5;
    H=(repmat(Hgage,1,LinkNum)./(Agage.^0.4)).*repmat(usarea',timesteps,1).^0.4;
    Btmax=max(B,[],1)';
    
    # Run our component
    nst.run_one_step

# -*- coding: utf-8 -*-
"""
Created on Thu Apr 16 11:25:10 2020

@author: jczuba
"""



#%% IMPORT

import warnings
#warnings.filterwarnings('ignore')

import os
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from landlab.components import FlowDirectorSteepest, NetworkSedimentTransporter
from landlab.data_record import DataRecord
from landlab.grid.network import NetworkModelGrid
from landlab.plot import graph
from landlab.io import read_shapefile

from landlab.plot import plot_network_and_parcels
#%matplotlib inline



#%% LOAD SHAPEFILE

#%% Methow -- Allison nodes
DATA_DIR = "C:/Users/jczuba/Documents/GitHub/landlab/landlab/data/io/shapefile/methow/"
file = os.path.join(DATA_DIR, "MethowSubBasin.shp")

points_shapefile = os.path.join(DATA_DIR, "MethowSubBasin_Nodes_4.shp")
grid = read_shapefile(
    file,
    points_shapefile=points_shapefile,
    node_fields=["usarea_km2", "Elev_m"],
    link_fields=["usarea_km2", "Length_m"],
    link_field_conversion={"usarea_km2": "drainage_area", "Slope":"channel_slope", "Length_m":"reach_length"},
    node_field_conversion={
        "usarea_km2": "drainage_area",
        "Elev_m": "topographic__elevation",
    },
    threshold=0.01,
    )

#%% Logan River
DATA_DIR = "./data/"
file = os.path.join(DATA_DIR, "a004_network.shp")

points_shapefile = os.path.join(DATA_DIR, "a004_nodes_att.shp")
grid = read_shapefile(
    file,
    points_shapefile=points_shapefile,
    node_fields=["usarea_km2", "elev_m"],
    link_fields=["usarea_km2", "Length_m"],
    link_field_conversion={"usarea_km2": "drainage_area", "Slope":"channel_slope", "Length_m":"reach_length"},
    node_field_conversion={
        "usarea_km2": "drainage_area",
        "elev_m": "topographic__elevation",
    },
    threshold=0.01,
    )

#%% 001 Red Butte
DATA_DIR = "C:/Users/jczuba/Documents/Projects/UtahWildfire/GIS/analysis/Watersheds/001/"

file = os.path.join(DATA_DIR, "a001_network.shp")

points_shapefile = os.path.join(DATA_DIR, "delete/a001_nodes_att.shp")
grid = read_shapefile(
    file,
    points_shapefile=points_shapefile,
    node_fields=["usarea_km2", "elev_m"],
    link_fields=["usarea_km2", "Length_m"],
    link_field_conversion={"usarea_km2": "drainage_area", "Slope":"channel_slope", "Length_m":"reach_length"},
    node_field_conversion={
        "usarea_km2": "drainage_area",
        "elev_m": "topographic__elevation",
    },
    threshold=0.01,
    )

#%%
grid.at_link.keys()
grid.at_node.keys()

graph.plot_graph(grid, at="node,link")

grid.number_of_links
grid.number_of_nodes

grid.at_node["bedrock__elevation"] = grid.at_node["topographic__elevation"].copy()



#%% FLOW

DATA_DIR = "C:/Users/jczuba/Documents/Python_local/network_sediment_transporter/data/"
file = os.path.join(DATA_DIR, "USGS 10172200 RED BUTTE CREEK AT FORT DOUGLAS, NEAR SLC, UT.xlsx")

data = pd.read_excel (
        file,
        sheet_name='Mean Daily Flow',
        header = 3)

#drop other columns if needed
#data.dropna(subset=['Stop Number','Lithology Category', 'Size (cm)'])

flow_at_gage_cfs = data['Discharge, ft3/s (mean)'].values

#convert Q in ft3/s to m3/s
flow_at_gage_cms = flow_at_gage_cfs * 0.3048 * 0.3048 * 0.3048 #%ft3/s to m3/s

#HELP HELP -- need to work on handling dates
#date = data['Date'].values

#plot hydrograph
#add date as x-axis
plt.plot(flow_at_gage_cms)
#add ylabel('Daily streamflow, m^3/s')

#HELP HELP -- identify and remove bad values
#Q(Q==0)=NaN;

#specify recurrence interval for which to calculate flow
return_period_years = 2 #years

# calculate recurrence interval flow
recurrence_interval_flow = calc_recurrence_interval_flow(
    flow_at_gage_cms, 
    return_period_years
    )

#%%

grid.at_link["channel_width"] = 1 * np.ones(grid.number_of_links) # m
grid.at_link["flow_depth"] = 0.5 * np.ones(grid.number_of_links) # m



#%% CREATE PARCELS

# element_id is the link on which the parcel begins. 
element_id = np.repeat(np.arange(grid.number_of_links), 1)
element_id = np.expand_dims(element_id, axis=1)

volume = 1*np.ones(np.shape(element_id))  # (m3)
active_layer = np.ones(np.shape(element_id)) # 1= active, 0 = inactive
density = 2650 * np.ones(np.size(element_id))  # (kg/m3)
abrasion_rate = 0 * np.ones(np.size(element_id)) # (mass loss /m)

# # Lognormal GSD
# medianD = 0.15 # m
# mu = np.log(medianD)
# sigma = np.log(2) #assume that D84 = sigma*D50
# np.random.seed(0)
# D = np.random.lognormal(
#     mu,
#     sigma,
#     np.shape(element_id)
# )  # (m) the diameter of grains in each parcel

D = 1.1 * np.ones(np.shape(element_id)) # m

time_arrival_in_link = np.random.rand(np.size(element_id), 1) 
#location_in_link = np.random.rand(np.size(element_id), 1) 
location_in_link = 0 * np.ones(np.shape(element_id))

lithology = ["quartzite"] * np.size(element_id)

variables = {
    "abrasion_rate": (["item_id"], abrasion_rate),
    "density": (["item_id"], density),
    "lithology": (["item_id"], lithology),
    "time_arrival_in_link": (["item_id", "time"], time_arrival_in_link),
    "active_layer": (["item_id", "time"], active_layer),
    "location_in_link": (["item_id", "time"], location_in_link),
    "D": (["item_id", "time"], D),
    "volume": (["item_id", "time"], volume)
}

items = {"grid_element": "link", "element_id": element_id}

_OUT_OF_NETWORK = NetworkModelGrid.BAD_INDEX - 1

parcels = DataRecord(
    grid,
    items=items,
    time=[0.0],
    data_vars=variables,
    dummy_elements={"link": [_OUT_OF_NETWORK]},
)

#%%
"""
But lets say you had a bumpy array called flow that was nlinks x n time steps. then you would just do. 

for i in range(flow.shape[1]):
   grid.at_link['flow_depth][:] = flow[:, i]
   nst.run_one_step(dt)

"""

#%% RUN NST

timesteps = 10 # total number of timesteps
dt = 60 * 60 * 24 *2 # length of timestep (seconds) 

fd = FlowDirectorSteepest(grid, "topographic__elevation")
fd.run_one_step()

nst = NetworkSedimentTransporter(    
    grid,
    parcels,
    fd,
    bed_porosity=0.3,
    g=9.81,
    fluid_density=1000,
    transport_method="WilcockCrowe",
)

for t in range(0, (timesteps * dt), dt):
    nst.run_one_step(dt)
    
    print("Model time: ", t/(60*60*24), "days passed")
    


#%% PLOT
        
fig = plot_network_and_parcels(grid, 
                               parcels,  
                               parcel_time_index=4, 
                               #parcel_filter=parcel_filter,
                               #link_attribute="sediment_total_volume", 
                               #network_norm=network_norm,
                               network_linewidth=4,
                               #network_cmap='bone_r',
                               parcel_alpha=1.0, 
                               parcel_color_attribute="D",
                               #parcel_color_norm=parcel_color_norm2, 
                               parcel_size_attribute="D",
                               parcel_size_min=5,
                               parcel_size_max=150,
                               #parcel_size_norm=parcel_size_norm,
                               parcel_size_attribute_title="D")

#%%
parcel_vol_on_grid = parcels.dataset["volume"].values
parcel_vol_on_grid[parcels.dataset["element_id"].values==-2]=0

#plt.figure(figsize=(8,6))
plt.plot(np.asarray(parcels.time_coordinates)/(60*60*24), 
         np.sum(parcel_vol_on_grid, axis=0),
         '-',
         linewidth=3, 
         alpha=0.5
        )

plt.ylabel('Total volume of parcels on grid $[m^3]$')
plt.xlabel('Time [days]')
plt.show() 

#%%

plt.loglog(parcels.dataset.D[:,-1],
         nst._distance_traveled_cumulative,
         '.'
        )
plt.xlabel('Parcel grain size (m)')
plt.ylabel('Cumulative parcel travel distance (m)')

#%%

plt.plot(grid.at_link["channel_slope"],
         nst.d_mean_active, 
         '.')
plt.xlabel('Channel slope (m/m)')
plt.ylabel('Mean grain size of active layer (m)')

#%%

plt.plot(grid.at_link["channel_slope"],
         nst.d_mean_active, 
         '.')
plt.xlabel('Channel slope (m/m)')
plt.ylabel('Mean grain size of active layer (m)')


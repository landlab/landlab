# -*- coding: utf-8 -*-
"""
This code outlines very a very basic use case for the NetworkSedimentTransporter
component, including adding new sediment mid-run

@authors: AMP
"""
import matplotlib.pyplot as plt
import numpy as np
import time as measuretime

from landlab.components import FlowDirectorSteepest, NetworkSedimentTransporter
from landlab.data_record import DataRecord
from landlab.grid.network import NetworkModelGrid

_OUT_OF_NETWORK = NetworkModelGrid.BAD_INDEX - 1

# %% Set the geometry using Network model grid (should be able to read in a shapefile here)

x_of_node = np.arange(0,1000000,10000)
y_of_node = np.zeros_like(x_of_node)

nodes_at_link = []
for i in range(np.size(y_of_node)-1): 
    nodes_at_link.append((i,i+1))

grid = NetworkModelGrid((y_of_node, x_of_node), nodes_at_link)
#plt.figure(0)
#graph.plot_graph(grid, at="node,link")

grid.at_node["topographic__elevation"] = 10000 - x_of_node*0.01

grid.at_node["bedrock__elevation"] = 10000 - x_of_node*0.01

area = grid.add_ones("cell_area_at_node", at="node")

# %% Set geometry for each link

grid.at_link["drainage_area"] = np.ones(grid.number_of_links) # m2
grid.at_link["channel_slope"] = np.ones(grid.number_of_links)
grid.at_link["reach_length"] = 10000*np.ones(grid.number_of_links)  # m
grid.at_link["channel_width"] = 1*np.ones(grid.number_of_links)


# %% initialize bed sediment 

timesteps = 100
time = [0.0]  

element_id = np.repeat(np.arange(99),200)
element_id = np.expand_dims(element_id, axis=1)

np.random.seed(0)

time_arrival_in_link = np.random.rand(np.size(element_id), 1) 
starting_link = np.squeeze(element_id)
volume = 5*np.ones(np.shape(element_id))  
lithology = ["quartzite"] * np.size(element_id)
active_layer = np.ones(np.shape(element_id))
density = 2650 * np.ones(np.size(element_id))  # (kg/m3)
location_in_link = np.random.rand(np.size(element_id), 1) 

# ------ ABRASION RATE: no abrasion
abrasion_rate = 0 * np.ones(np.size(element_id))

# ------ GRAIN SIZE:Lognormal GSD ------
medianD = 0.085
mu = np.log(medianD)
sigma = np.log(1.5) #assume that D84 = 2.3*D50
D = np.random.lognormal(mu,sigma,np.shape(element_id))  # (m) the diameter of grains in each parcel

items = {"grid_element": "link", "element_id": element_id}

variables = {
    "starting_link": (["item_id"], starting_link),
    "abrasion_rate": (["item_id"], abrasion_rate),
    "density": (["item_id"], density),
    "lithology": (["item_id"], lithology),
    "time_arrival_in_link": (["item_id", "time"], time_arrival_in_link),
    "active_layer": (["item_id", "time"], active_layer),
    "location_in_link": (["item_id", "time"], location_in_link),
    "D": (["item_id", "time"], D),
    "volume": (["item_id", "time"], volume),
}

parcels = DataRecord(
    grid,
    items=items,
    time=time,
    data_vars=variables,
    dummy_elements={"link": [_OUT_OF_NETWORK]},
)

# %% Flow parameters 

dt = 60 * 60 * 24 *2 # (seconds) timestep
flow_depth =  2 * np.ones([timesteps + 1,grid.number_of_links])

# %% Instantiate component(s)

fd = FlowDirectorSteepest(grid, "topographic__elevation")
fd.run_one_step()

nst = NetworkSedimentTransporter(
        
    grid,
    parcels,
    fd,
    flow_depth,
    bed_porosity=0.3,
    g=9.81,
    fluid_density=1000,
    transport_method="WilcockCrowe",
)

# %% Run the component(s)
pulsetime = 1
num_pulse_parcels = 2
Pulse_cum_transport = []

for t in range(0, (timesteps * dt), dt):
    start_time = measuretime.time()
    #print("timestep ", [t], "started")
    # Run our component
    nst.run_one_step(dt)
    
    # RECYCLE sediment: what left the network gets added back in at top. 
    parcels.dataset.element_id.values[parcels.dataset.element_id.values == _OUT_OF_NETWORK] = 0
    
    if t == dt*pulsetime:
        print("Now let's add new parcels!")
        
        newpar_element_id = np.zeros(num_pulse_parcels,dtype=int)
        newpar_element_id = np.expand_dims(newpar_element_id, axis=1)
        
        new_starting_link = np.squeeze(newpar_element_id)
        
        np.random.seed(0)
        
        new_time_arrival_in_link = nst._time* np.ones(
            np.shape(newpar_element_id)) 
        
        new_volume = 0.5*np.ones(np.shape(newpar_element_id))  # (m3) the volume of each parcel
        
        new_lithology = ["pulse_material"] * np.size(
            newpar_element_id
        )  # a lithology descriptor for each parcel
        
        new_active_layer = np.ones(
            np.shape(newpar_element_id)
        )  # 1 = active/surface layer; 0 = subsurface layer
        
        new_density = 2650 * np.ones(np.size(newpar_element_id))  # (kg/m3)
        
        new_location_in_link = np.random.rand(
            np.size(newpar_element_id), 1
        )

        new_abrasion_rate = 0 * np.ones(np.size(newpar_element_id))
        
        new_D = 0.03 * np.ones(np.shape(newpar_element_id))

        newpar_grid_elements = np.array(
            np.empty(
                (np.shape(newpar_element_id)), dtype=object
            )
        ) # BUG: should be able to pass ["link"], but datarecord fills it into an incorrect array shape-- the length of parcels (NOT new parcels)
        newpar_grid_elements.fill("link")
        
        new_parcels = {"grid_element": newpar_grid_elements,
                 "element_id": newpar_element_id}

        new_variables = {
            "starting_link": (["item_id"], new_starting_link),
            "abrasion_rate": (["item_id"], new_abrasion_rate),
            "density": (["item_id"], new_density),
            "lithology": (["item_id"], new_lithology),
            "time_arrival_in_link": (["item_id", "time"], new_time_arrival_in_link),
            "active_layer": (["item_id", "time"], new_active_layer),
            "location_in_link": (["item_id", "time"], new_location_in_link),
            "D": (["item_id", "time"], new_D),
            "volume": (["item_id", "time"], new_volume),
        }
        
        parcels.add_item(
                time=[nst._time],
                new_item = new_parcels,
                new_item_spec = new_variables
        )
        
        Pulse_cum_transport = np.zeros(
                [len(nst._distance_traveled_cumulative[-num_pulse_parcels:]),
                 1
                 ]
                )
    
    if t > dt*pulsetime :
        Pulse_cum_transport = np.append(
                Pulse_cum_transport,
                np.expand_dims(
                        nst._distance_traveled_cumulative[-num_pulse_parcels:],
                        axis = 1),
                axis = 1
                )
    
    print("Model time: ", t/(60*60*24), "days passed")
    print('Elapsed:', (measuretime.time() - start_time)/60 ,' minutes')
    print()

plt.show()

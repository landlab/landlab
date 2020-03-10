# -*- coding: utf-8 -*-
"""
This code outlines very a very basic use case for the NetworkSedimentTransporter
component.

Created on Sun May 20 15:54:03 2018

@authors: Jon Czuba, Allison Pfeiffer, Katy Barnhart
"""
import matplotlib.pyplot as plt
import numpy as np

# from landlab.components import NetworkSedimentTransporter
from landlab.components import FlowDirectorSteepest, NetworkSedimentTransporter
from landlab.data_record import DataRecord
from landlab.grid.network import NetworkModelGrid
from landlab.plot import graph
from landlab.plot.network_sediment_transporter import (
    plot_network_links,
    plot_network_parcels,
    plot_pathway_values,
)

# from landlab.plot.network_sediment_transporter import *  # Note-- this is an example. it loads plotting scripts that don't exist yet.

_OUT_OF_NETWORK = NetworkModelGrid.BAD_INDEX - 1

# %% Set the geometry using Network model grid (should be able to read in a shapefile here)

y_of_node = (0, 100, 200, 200, 300, 400, 400, 125)
x_of_node = (0, 0, 100, -50, -100, 50, -150, -100)

nodes_at_link = ((1, 0), (2, 1), (1, 7), (3, 1), (3, 4), (4, 5), (4, 6))

grid = NetworkModelGrid((y_of_node, x_of_node), nodes_at_link)
plt.figure(0)
graph.plot_graph(grid, at="node,link")

grid.at_node["topographic__elevation"] = [0.0, 0.1, 0.3, 0.3, 0.45, 0.45, 0.7, 0.8]
grid.at_node["bedrock__elevation"] = [0.0, 0.1, 0.3, 0.3, 0.45, 0.45, 0.7, 0.8]

area = grid.add_ones("cell_area_at_node", at="node")

# %% Set geometry for each link

# Ultimately, map between flow accumulator and shapefile reader info...
# map_upstream_node_to_link

grid.at_link["drainage_area"] = [100e6, 10e6, 70e6, 20e6, 70e6, 30e6, 40e6]  # m2
grid.at_link["channel_slope"] = [0.01, 0.02, 0.01, 0.02, 0.02, 0.03, 0.03]
grid.at_link["reach_length"] = [
    10000.0,
    10000.0,
    10000.0,
    10000.0,
    10000.0,
    10000.0,
    10000.0,
]  # m
## Why are we not using grid.length_of_link? Is it because this is the full
# squiggly path rather than the node-to-node.


grid.at_link["channel_width"] = 15 * np.ones(
    np.size(grid.at_link["drainage_area"])
)  # m REPLACE with something hydraulically meaningful
grid.at_link["channel_width"][3] = 10
# modify elevations so they are consistent with adjusted slopes

## Basic parameters

g = 9.81  # m/s2
rho = 1000  # kg/m3

bed_porosity = 0.3  # porosity of the bed material

# %% initialize bed sediment (may become its own component)

# NOTE: inputs to DataRecord need to have the same shape as the time/item inputs
# So, if a parcel attribute is being tracked in time, it needs to have
# np.shape = (n,1), and if it isn't tracked in time, it needs to have
# np.shape = (n,).

# Ultimately,
# parcels = SedimentParcels(grid,initialization_info_including_future_forcing)


timesteps = 8

element_id = np.repeat(np.array([0, 1, 2, 3, 4, 5, 6], dtype=int), 100)

#
# element_id = np.array(
#    [0, 0, 1, 1, 1,
#     5, 2, 2, 3, 3, 4, 4, 5, 5, 5, 5, 6, 6, 2, 3, 4, 4, 4, 3, 4, 5],
#    dtype=int,
# ) # current link for each parcel

element_id = np.expand_dims(element_id, axis=1)

starting_link = np.squeeze(element_id)  # starting link for each parcel

np.random.seed(0)

time_arrival_in_link = np.random.rand(
    np.size(element_id), 1
)  # time of arrival in each link -- larger numbers are younger
volume = np.ones(np.shape(element_id))  # (m3) the volume of each parcel
D = 0.05 * np.ones(np.shape(element_id))  # (m) the diameter of grains in each parcel
lithology = ["quartzite"] * np.size(
    element_id
)  # a lithology descriptor for each parcel
abrasion_rate = 0.0001 * np.ones(
    np.size(element_id)
)  # 0 = no abrasion; abrasion rates are positive mass loss coefficients (mass loss / METER)
active_layer = np.ones(
    np.shape(element_id)
)  # 1 = active/surface layer; 0 = subsurface layer

density = 2650 * np.ones(np.size(element_id))  # (kg/m3)

location_in_link = np.zeros(
    np.shape(element_id)
)  # [0 1], 0 is upstream end of link, 1 is downstream end

D[0] = 0.075
D[5] = 0.0001  # make one of them sand

volume[2] = 0.3

time = [0.0]  # probably not the sensible way to do this...

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


# Add parcels in at a given time --> attribute in the item collection

# %% Flow parameters (this will happen in the sq.run_one_step and sc.run_one_step)

# Made up hydraulic geometry

Qgage = 1000.0  # (m3/s)
dt = 60 * 60 * 24  # (seconds) daily timestep

Bgage = 30.906 * Qgage ** 0.1215
# (m)
Hgage = 1.703 * Qgage ** 0.3447
# (m)
Agage = 4.5895e9
# (m2)

channel_width = (np.tile(Bgage, (grid.number_of_links)) / (Agage ** 0.5)) * np.tile(
    grid.at_link["drainage_area"], (timesteps, 1)
) ** 0.5

flow_depth = (np.tile(Hgage, (grid.number_of_links)) / (Agage ** 0.4)) * np.tile(
    grid.at_link["drainage_area"], (timesteps + 1, 1)
) ** 0.4

Btmax = np.amax(channel_width, axis=0)  # CURRENTLY UNUSED

# %% Instantiate component(s)
# dis = ExteralDischargeSetter(grid, filename, model='dhsvm')

# sq = SyntheticDischargeMaker(discharge,drainage_area) # OR read DHSVM.
# define a surface_water_discharge for each link and each timestep

# sc = SyntheticChannelGeomMaker(hydraulic_geometry_scaling_rules,discharge)
#

fd = FlowDirectorSteepest(grid, "topographic__elevation")
fd.run_one_step()

nst = NetworkSedimentTransporter(
    grid,
    parcels,
    fd,
    flow_depth,
    bed_porosity,
    g=9.81,
    fluid_density=1000,
    transport_method="WilcockCrowe",
)

# %% Run the component(s)

for t in range(0, (timesteps * dt), dt):
    # print("timestep ", [t], "started")
    print("Model time: ", t / (60 * 60 * 24), "days passed")
    # move any sediment additions from forcing Item collector to bed item collector

    # sq.run_one_step
    #   will assign discharge for each reach

    # sc.run_one_step
    #   will assign flow depth for each reach (for this timestep)

    # Run our component
    nst.run_one_step(dt)


#%%
# from landlab.plot.network import plot_pathway
# plot_pathway(nst, time=3, link=4)

# %% A few plot outputs, just to get started.

plt.figure(1)
plt.plot(parcels.time_coordinates, parcels.dataset.location_in_link.values[0, :], ".")
plt.plot(parcels.time_coordinates, parcels.dataset.location_in_link.values[14, :], ".")
plt.plot(parcels.time_coordinates, parcels.dataset.location_in_link.values[16, :], ".")
plt.plot(parcels.time_coordinates, parcels.dataset.location_in_link.values[17, :], ".")

plt.title("Tracking link location for a single parcel")
plt.xlabel("time")
plt.ylabel("location in link")


plt.figure(2)
plt.plot(parcels.time_coordinates, parcels.dataset.element_id.values[0, :], ".")
plt.plot(parcels.time_coordinates, parcels.dataset.element_id.values[14, :], ".")
plt.plot(parcels.time_coordinates, parcels.dataset.element_id.values[16, :], ".")
plt.plot(parcels.time_coordinates, parcels.dataset.element_id.values[17, :], ".")

plt.title("Tracking link location for a single parcel")
plt.xlabel("time")
plt.ylabel("link")


plt.figure(3)
plt.plot(
    parcels.time_coordinates, np.sum(parcels.dataset["volume"].values, axis=0), "."
)
plt.title("Silly example: total volume, all parcels through time")
plt.xlabel("time")
plt.ylabel("total volume of parcels")

# %% Plotting routines to step through
# %%
# %% 1. Run "import_shapefile.py" ideally, we should be reading a shapefile at the outset of NST_setup to gather our network properties from the shapefile attribute table. This would also include loading in, early on, the squiggly line network. For now this needs to be done here before plotting because it has not been done yet.

DATA_DIR = "./data/"

file = os.path.join(DATA_DIR, "Test_Network.shp")
grid = read_shapefile(file)

# %% 2. Plot an attribute on the delineated network

# determine link attribute to plot

# This is really a time index to plot, so likely change the name
parcel_time = 5

# NEXT NEXT: Put these next few lines into a function: composite_parcel_attributes_at_links
# Also make the inverse: map_link_attributes_to_parcels
# Initialize array of a composite link attribute
lnkvol = np.empty(grid.at_link.size)
lnkvol[:] = np.nan
# Composite link attributes
for i in range(grid.at_link.size):
    lnkvol[i] = np.sum(
        parcels.dataset.volume[
            (parcels.dataset.element_id[:, parcel_time]) == i, parcel_time
        ]
    )
link_attribute = lnkvol
# link_attribute = grid.at_link["GridID"]

# plot_network_links(nmg, nmg.at_link["channel_width"])
# plot_network_links(nmg, an_array_I_made_from_fancy_calculations)
plot_network_links(grid, link_attribute)

# %% 3. Plot parcels as dots on the delineated network

# NOTE: We could also easily plot a single parcel at all times. Currently, this code is just set up to plot all parcels at a single time.

# This is really a time index to plot, so likely change the name
parcel_time = 8

# Determine color and size of parcels based on a parcel attribute
parcel_color = parcels.dataset.element_id[:, parcel_time]
parcel_size = 15

# Determine color and size of parcels based on a link attribute
# Need to map link attributes back to parcels based on element id
# NEXT NEXT: based on the code above it will be worthwhile to write a function that maps link attributes to parcels and another function that compiles link attributes from parcels. Then we can just call that function above and here.
# parcel_color = parcels.dataset.element_id[:,parcel_time]
# parcel_size = 15

plot_network_parcels(grid, parcels, parcel_time, parcel_color, parcel_size)

# %% 4. Plot link attributes along a pathway through the network.

# determine a specific link on the network
link_number = 4
# determine an attribute on the grid to plot
link_attribute = "GridID"

plot_pathway_values(link_number, link_attribute, fd, grid)


# %% 5. Develop GUI to more quickly interface with calling these plotting functions.

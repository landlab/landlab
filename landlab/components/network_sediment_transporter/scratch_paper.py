# -*- coding: utf-8 -*-
"""
Created on Thu Sep 12 12:06:29 2019

Temporary file for use while writing tests


@author: pfeif
"""
import numpy as np
import matplotlib.pyplot as plt
from landlab.components import FlowDirectorSteepest, NetworkSedimentTransporter
from landlab.grid.network import NetworkModelGrid
from landlab import BAD_INDEX_VALUE
from landlab.data_record import DataRecord
from landlab.plot import graph
_OUT_OF_NETWORK = BAD_INDEX_VALUE - 1

#Create a network model grid to represent the channel network
y_of_node = (0, 0, 0, 0)
x_of_node = (0, 100, 200, 300)
nodes_at_link = ((0,1), (1,2), (2,3))

nmg = NetworkModelGrid((y_of_node, x_of_node), nodes_at_link)

plt.figure(0)
graph.plot_graph(nmg, at="node,link")

# %%
# Add variables to the NetworkModelGrid
nmg.at_node["topographic__elevation"] = [3., 2., 1., 0.] # m
nmg.at_node["bedrock__elevation"] = [3., 2., 1., 0.] # m
nmg.at_link["drainage_area"] = [10e6, 10e6, 10e6]  # m2
nmg.at_link["channel_slope"] = [0.001, 0.001, 0.001]
nmg.at_link["link_length"] = [100, 100, 100]  # m
nmg.at_link["channel_width"] = (15 * np.ones(np.size(nmg.at_link["drainage_area"])))

flow_director = FlowDirectorSteepest(nmg)
flow_director.run_one_step()

timesteps = 10

example_flow_depth = (
    np.tile(2, (nmg.number_of_links))
) * np.tile(1, (timesteps + 1, 1)) # 2 meter flow depth

time = [0.0]

# Set up sediment parcels DataRecord
items = {"grid_element": "link",
         "element_id": np.array([[0]])}

variables = {
    "starting_link": (["item_id"], np.array([0])),
    "abrasion_rate": (["item_id"], np.array([0])),
    "density": (["item_id"], np.array([2650])),
    "time_arrival_in_link": (["item_id", "time"], np.array([[0]])),
    "active_layer": (["item_id", "time"], np.array([[1]])),
    "location_in_link": (["item_id", "time"], np.array([[0]])),
    "D": (["item_id", "time"], np.array([[0.05]])),
    "volume": (["item_id", "time"], np.array([[1]])),
}

one_parcel = DataRecord(
    nmg,
    items=items,
    time=time,
    data_vars=variables,
    dummy_elements={"link": [_OUT_OF_NETWORK]},
)
# Instantiate the model run
nst = NetworkSedimentTransporter(
        nmg,
        one_parcel,
        flow_director,
        example_flow_depth,
        bed_porosity=0.03,
        g=9.81,
        fluid_density=1000,
        channel_width="channel_width",
        transport_method="WilcockCrowe",
    )

dt = 60  # (seconds) 1 min timestep

# Run the model
for t in range(0, (timesteps * dt), dt):
        nst.run_one_step(dt)



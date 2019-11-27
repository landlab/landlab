# -*- coding: utf-8 -*-
"""
Created on Thu Sep 12 12:06:29 2019

Temporary file for use while writing tests


@author: pfeif
"""

import matplotlib.pyplot as plt
import numpy as np

# from landlab.components import NetworkSedimentTransporter
from landlab import BAD_INDEX_VALUE
from landlab.components import FlowDirectorSteepest, NetworkSedimentTransporter
from landlab.data_record import DataRecord
from landlab.grid.network import NetworkModelGrid
from landlab.plot import graph
from numpy.testing import assert_array_almost_equal, assert_array_equal

_OUT_OF_NETWORK = BAD_INDEX_VALUE - 1

# %%  


#
y_of_node = (0, 0, 0, 0)
x_of_node = (0, 100, 200, 300)
nodes_at_link = ((0,1), (1,2), (2,3))

nmg_constant_slope = NetworkModelGrid((y_of_node, x_of_node), nodes_at_link)

plt.figure(0)
graph.plot_graph(nmg_constant_slope, at="node,link")

# add variables to nmg
nmg_constant_slope.at_node["topographic__elevation"] = [3., 2., 1., 0.]
nmg_constant_slope.at_node["bedrock__elevation"] = [3., 2., 1., 0.]
area = nmg_constant_slope.add_ones("cell_area_at_node", at="node")
nmg_constant_slope.at_link["drainage_area"] = [10e6, 10e6, 10e6]  # m2
nmg_constant_slope.at_link["channel_slope"] = [0.001, 0.001, 0.001]
nmg_constant_slope.at_link["link_length"] = [100, 100, 100]  # m

nmg_constant_slope.at_link["channel_width"] = 15 * np.ones(np.size(nmg_constant_slope.at_link["drainage_area"]))

flow_director = FlowDirectorSteepest(nmg_constant_slope)
flow_director.run_one_step()

timesteps = 11
#timesteps = 33 # playing to demonstrate issue...

example_flow_depth = (
    np.tile(2, (nmg_constant_slope.number_of_links))
) * np.tile(1, (timesteps + 1, 1))
# 2 meter flow depth

#example_flow_depth = example_flow_depth*0.5 

time = [0.0]  # probably not the sensible way to do this...

element_id = np.zeros(100, dtype=int)
element_id = np.expand_dims(element_id, axis=1)

items = {"grid_element": "link",
         "element_id": element_id}

abrasion_rate = np.zeros(np.size(element_id))
initial_volume = np.ones(np.shape(element_id))
time_arrival_in_link = np.arange(0,0.1,0.001)
time_arrival_in_link = np.expand_dims(time_arrival_in_link,axis=1)

variables = {
    "starting_link": (["item_id"], np.squeeze(element_id)),
    "abrasion_rate": (["item_id"], abrasion_rate),
    "density": (["item_id"], 2650*np.ones(np.size(element_id))),
    "time_arrival_in_link": (["item_id", "time"], time_arrival_in_link),
    "active_layer": (["item_id", "time"],initial_volume),
    "location_in_link": (["item_id", "time"], element_id),
    "D": (["item_id", "time"], initial_volume*0.05),
    "volume": (["item_id", "time"], initial_volume),
}

hundred_boring_parcels = DataRecord(
    nmg_constant_slope,
    items=items,
    time=time,
    data_vars=variables,
    dummy_elements={"link": [_OUT_OF_NETWORK]},
)

#example_flow_depth = example_flow_depth*5# outrageously high transport rate

nst = NetworkSedimentTransporter(
        nmg_constant_slope,
        hundred_boring_parcels,
        flow_director,
        example_flow_depth,
        bed_porosity=0.03,
        g=9.81,
        fluid_density=1000,
        channel_width="channel_width",
        transport_method="WilcockCrowe",
    )

dt = 60 *60 # (seconds) 1 hour timestep

distance_traveled = np.arange(0.,timesteps)
pvelocity = np.arange(0.,timesteps)
active_layer_thickness_array = np.arange(0.,timesteps)


for t in range(0, (timesteps * dt), dt):
        nst.run_one_step(dt)
        distance_traveled[np.int(t/dt)]= nst._distance_traveled_cumulative[-1]
        print('distance_traveled_cumulative', nst._distance_traveled_cumulative)
        active_layer_thickness_array[np.int(t/dt)]=nst.active_layer_thickness_array[-1]

#print("Cumulative distance traveled", 
#      nst._distance_traveled_cumulative)

#assert_array_almost_equal(x,y)

plt.plot(hundred_boring_parcels.time_coordinates, 
         hundred_boring_parcels.dataset.location_in_link.values[-1, :], ".")

#plt.plot(active_layer_thickness_array,'.')
plt.plot(hundred_boring_parcels.dataset.time_arrival_in_link.values[-1,:],'.')
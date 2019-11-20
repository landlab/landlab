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

items = {"grid_element": "link",
         "element_id": np.array([[0]])}

initial_volume = np.array([[1]])
abrasion_rate = np.array([0])

variables = {
    "starting_link": (["item_id"], np.array([0])),
    "abrasion_rate": (["item_id"], abrasion_rate),
    "density": (["item_id"], np.array([2650])),
    "time_arrival_in_link": (["item_id", "time"], np.array([[0.71518937]])),
    "active_layer": (["item_id", "time"], np.array([[1]])),
    "location_in_link": (["item_id", "time"], np.array([[0]])),
    "D": (["item_id", "time"], np.array([[0.05]])),
    "volume": (["item_id", "time"], initial_volume),
}

one_parcel = DataRecord(
    nmg_constant_slope,
    items=items,
    time=time,
    data_vars=variables,
    dummy_elements={"link": [_OUT_OF_NETWORK]},
)

#example_flow_depth = example_flow_depth*5# outrageously high transport rate

nst = NetworkSedimentTransporter(
        nmg_constant_slope,
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

distance_traveled = np.arange(0.,timesteps)
pvelocity = np.arange(0.,timesteps)
active_layer_thickness_array = np.arange(0.,timesteps)
#distance_traveled = np.arange(0.,timesteps)

for t in range(0, (timesteps * dt), dt):
        nst.run_one_step(dt)
        distance_traveled[np.int(t/dt)]= nst._distance_traveled_cumulative
        active_layer_thickness_array[np.int(t/dt)]=nst.active_layer_thickness_array
# COMPARE TO: 
        # 2 meter flow depth
        # 0.01 slope 
        # link length = 100 m
        # 0.05 m grain size parcel
        # 1 m3 parcel volume
        # timestep = 1 min
        # notes in notebook. transport distance in 1 min =21.619985271513052
        # total transport distance = 237.81983798664356 m 
        # final location in link after 11 timesteps: 0.3781983798664356
 
print("Parcel location in link", 
      nst._parcels.dataset.location_in_link[0,-1])

#print("Cumulative distance traveled", 
#      nst._distance_traveled_cumulative)

#assert_array_almost_equal(x,y)

plt.plot(one_parcel.time_coordinates, 
         one_parcel.dataset.location_in_link.values[0, :], ".")

plt.plot(active_layer_thickness_array,'.')

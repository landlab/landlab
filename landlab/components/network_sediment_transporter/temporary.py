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



y_of_node = (0, 100, 200, 200, 300, 400, 400, 125)
x_of_node = (0, 0, 100, -50, -100, 50, -150, -100)
nodes_at_link = ((1, 0), (2, 1), (1, 7), (3, 1), (3, 4), (4, 5), (4, 6))

example_nmg = NetworkModelGrid((y_of_node, x_of_node), nodes_at_link)

# add variables to example_nmg
example_nmg.at_node["topographic__elevation"] = [0.0, 0.1, 0.3, 0.2, 0.35, 0.45, 0.5, 0.6]
example_nmg.at_node["bedrock__elevation"] = [0.0, 0.1, 0.3, 0.2, 0.35, 0.45, 0.5, 0.6]
area = example_nmg.add_ones("cell_area_at_node", at="node")
example_nmg.at_link["drainage_area"] = [100e6, 10e6, 70e6, 20e6, 70e6, 30e6, 40e6]  # m2
example_nmg.at_link["channel_slope"] = [0.01, 0.02, 0.01, 0.02, 0.02, 0.03, 0.03]
example_nmg.at_link["link_length"] = [10000, 10000, 10000, 10000, 10000, 10000, 10000]  # m

example_nmg.at_link["channel_width"] = 15 * np.ones(np.size(example_nmg.at_link["drainage_area"]))


element_id = np.array(
    [0, 0, 1, 1, 2, 2, 3, 3, 4, 4, 5, 5, 6, 6, 6],
    dtype=int,
)  # current link for each parcel

element_id = np.expand_dims(element_id, axis=1)
starting_link = np.squeeze(element_id)  # starting link for each parcel

np.random.seed(0)

time_arrival_in_link = np.random.rand(
    np.size(element_id), 1
)  # time of arrival in each link -- larger numbers are younger
volume = np.ones(np.shape(element_id))  # (m3) the volume of each parcel
D = 0.05 * np.ones(
    np.shape(element_id)
)  # (m) the diameter of grains in each parcel
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
    "time_arrival_in_link": (["item_id", "time"], time_arrival_in_link),
    "active_layer": (["item_id", "time"], active_layer),
    "location_in_link": (["item_id", "time"], location_in_link),
    "D": (["item_id", "time"], D),
    "volume": (["item_id", "time"], volume),
}

example_parcels = DataRecord(
    example_nmg,
    items=items,
    time=time,
    data_vars=variables,
    dummy_elements={"link": [_OUT_OF_NETWORK]},
)

example_flow_director = FlowDirectorSteepest(example_nmg)
example_flow_director.run_one_step()

timesteps = 30

Qgage = 80000.0  # (m3/s)
dt = 60 * 60 * 24  # (seconds) daily timestep

Hgage = 1.703 * Qgage ** 0.3447
# (m)
Agage = 4.5895e9
# (m2)

example_flow_depth = (
    np.tile(Hgage, (example_nmg.number_of_links)) / (Agage ** 0.4)
) * np.tile(example_nmg.at_link["drainage_area"], (timesteps + 1, 1)) ** 0.4


# %% Ok, now to the tests. 
### Make a single parcel

time = [0.0]  # probably not the sensible way to do this...

items = {"grid_element": "link",
         "element_id": np.array([[6],[6]])}

initial_volume = np.array([[1],[1]])
abrasion_rate = np.array([0.0001,0])

variables = {
    "starting_link": (["item_id"], np.array([6,6])),
    "abrasion_rate": (["item_id"], abrasion_rate),
    "density": (["item_id"], np.array([2650,90000])),
    "time_arrival_in_link": (["item_id", "time"], np.array([[0.5],[0]])),
    "active_layer": (["item_id", "time"], np.array([[1],[1]])),
    "location_in_link": (["item_id", "time"], np.array([[0],[0]])),
    "D": (["item_id", "time"], np.array([[0.05],[0.05]])),
    "volume": (["item_id", "time"], initial_volume),
}

two_parcels = DataRecord(
    example_nmg,
    items=items,
    time=time,
    data_vars=variables,
    dummy_elements={"link": [_OUT_OF_NETWORK]},
)
    
timesteps = 8

example_flow_depth = example_flow_depth*5 # outrageously high transport rate

nst = NetworkSedimentTransporter(
        example_nmg,
        two_parcels,
        example_flow_director,
        example_flow_depth,
        bed_porosity=0.03,
        g=9.81,
        fluid_density=1000,
        channel_width="channel_width",
        transport_method="WilcockCrowe",
    )

dt = 60 * 60 * 24  # (seconds) daily timestep

for t in range(0, (timesteps * dt), dt):
    nst.run_one_step(dt)
    print("Successfully completed a timestep")
    # Need to define original_node_elev after a timestep has passed.
    if t/(60*60*24) == 1:
        original_parcel_vol = example_parcels.dataset.volume[:,0][0]
        

# Parcel volume should decrease according to abrasion rate
#distance_traveled = XXXX should calculate this properly... add up distance traveled in parcels as we go
#
#volume_after_transport = initial_volume * np.exp(distance_traveled*-abrasion_rate)
#assert_array_equal(
#        volume_after_transport,
#        two_parcels.dataset.volume[0:2,-1]
#        )




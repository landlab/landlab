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
from landlab import BAD_INDEX_VALUE
from landlab.components import FlowDirectorSteepest, NetworkSedimentTransporter
from landlab.data_record import DataRecord
from landlab.grid.network import NetworkModelGrid
from landlab.plot import graph

_OUT_OF_NETWORK = BAD_INDEX_VALUE - 1

# %% Set the geometry using Network model grid (should be able to read in a shapefile here)
def test_example(
    example_nmg, example_parcels, example_flow_director, example_flow_depth
):

    nst = NetworkSedimentTransporter(
        example_nmg,
        example_parcels,
        example_flow_director,
        example_flow_depth,
        active_layer_thickness=0.5,
        bed_porosity=0.03,
        g=9.81,
        fluid_density=1000,
        channel_width="channel_width",
        transport_method="WilcockCrowe",
    )

    # %% Run the component(s)
    dt = 60 * 60 * 24  # (seconds) daily timestep

    timesteps = 30
    for t in range(0, (timesteps * dt), dt):
        print("timestep ", [t], "started")
        # move any sediment additions from forcing Item collector to bed item collector

        # sq.run_one_step
        #   will assign discharge for each reach

        # sc.run_one_step
        #   will assign flow depth for each reach (for this timestep)

        # Run our component
        nst.run_one_step(dt)

    # %% A few plot outputs, just to get started.


plt.title("Silly example: total volume, all parcels through time")
plt.xlabel("time")
plt.ylabel("total volume of parcels")

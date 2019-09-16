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

def test_example(
    example_nmg, example_parcels, example_flow_director, example_flow_depth
):

    nst = NetworkSedimentTransporter(
        example_nmg,
        example_parcels,
        example_flow_director,
        example_flow_depth,
        bed_porosity=0.3,
        g=9.81,
        fluid_density=1000,
        channel_width="channel_width",
        transport_method="WilcockCrowe",
    )

    # %% Run the component(s)
    dt = 60 * 60 * 24  # (seconds) daily timestep

    timesteps = 30
    for t in range(0, (timesteps * dt), dt):
        nst.run_one_step(dt)

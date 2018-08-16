#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Feb 27 16:25:11 2018

@author: barnhark
"""
import numpy as np

from landlab import RasterModelGrid
from landlab.components import FlowAccumulator, FastscapeEroder, LinearDiffuser, DepressionFinderAndRouter
from landlab.plot import analyze_channel_network_and_plot

from nose.tools import assert_raises

def test_assertion_error():
    """Test that the correct assertion error will be raised."""

    mg = RasterModelGrid(10, 10)
    z = mg.add_zeros('topographic__elevation', at='node')
    z += 200 + mg.x_of_node + mg.y_of_node + np.random.randn(mg.size('node'))

    mg.set_closed_boundaries_at_grid_edges(bottom_is_closed=True, left_is_closed=True, right_is_closed=True, top_is_closed=True)
    mg.set_watershed_boundary_condition_outlet_id(0, z, -9999)
    fa = FlowAccumulator(mg, flow_director='D8', depression_finder=DepressionFinderAndRouter)
    sp = FastscapeEroder(mg, K_sp=.0001, m_sp=.5, n_sp=1)
    ld = LinearDiffuser(mg, linear_diffusivity=0.0001)

    dt = 100
    for i in range(200):
        fa.run_one_step()
        flooded = np.where(fa.depression_finder.flood_status==3)[0]
        sp.run_one_step(dt=dt,  flooded_nodes=flooded)
        ld.run_one_step(dt=dt)
        mg.at_node['topographic__elevation'][0] -= 0.001 # Uplift


    assert_raises(AssertionError,
                  analyze_channel_network_and_plot,
                  mg,
                  threshold = 100,
                  starting_nodes = [0],
                  number_of_channels=2)

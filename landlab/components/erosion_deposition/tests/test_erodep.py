#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Thu Jul 27 14:23:25 2017

@author: gtucker
"""

from landlab import RasterModelGrid, HexModelGrid
from landlab.components import ErosionDeposition, FlowAccumulator
import numpy as np
from numpy.testing import assert_equal
from nose.tools import assert_raises

def test_bad_solver_name():
    """
    Test that any solver name besides 'basic' and 'adaptive' raises an error.
    """

    #set up a 5x5 grid with one open outlet node and low initial elevations.
    nr = 5
    nc = 5
    mg = RasterModelGrid((nr, nc), 10.0)

    mg.add_zeros('node', 'topographic__elevation')

    mg['node']['topographic__elevation'] += mg.node_y / 10000 \
        + mg.node_x / 10000 \
        + np.random.rand(len(mg.node_y)) / 10000
    mg.set_closed_boundaries_at_grid_edges(bottom_is_closed=True,
                                           left_is_closed=True,
                                           right_is_closed=True,
                                           top_is_closed=True)
    mg.set_watershed_boundary_condition_outlet_id(0,
                                                  mg['node']['topographic__elevation'],
                                                  -9999.)

    # Create a D8 flow handler
    fa = FlowAccumulator(mg, flow_director='D8')

    #try to instantiate ErodionDeposition using a wrong solver name
    assert_raises(ValueError, ErosionDeposition, mg, K=0.01,
                         phi=0.0, v_s=0.001, m_sp=0.5, n_sp=1.0, 
                         sp_crit=0, F_f=0.0, solver='something_else')


def test_steady_state_with_basic_solver_option():
    pass


def test_can_run_with_hex():
    """Test that model can run with hex model grid."""

    # Set up a 5x5 grid with open boundaries and low initial elevations.
    mg = HexModelGrid(7, 7)
    z = mg.add_zeros('node', 'topographic__elevation')
    z[:] = 0.01 * mg.x_of_node

    # Create a D8 flow handler
    fa = FlowAccumulator(mg, flow_director='FlowDirectorSteepest')

    # Parameter values for test 1
    K = 0.001
    vs = 0.0001
    U = 0.001
    dt = 10.0

    # Create the ErosionDeposition component...
    ed = ErosionDeposition(mg, K=K, phi=0.0, v_s=vs, m_sp=0.5, n_sp=1.0,
                           solver='adaptive')

    # ... and run it to steady state.
    for i in range(2000):
        fa.run_one_step()
        ed.run_one_step(dt=dt)
        z[mg.core_nodes] += U * dt

    # Test the results
    s = mg.at_node['topographic__steepest_slope']
    sa_factor = (1.0 + vs) * U / K
    a18 = mg.at_node['drainage_area'][18]
    a28 = mg.at_node['drainage_area'][28]
    s = mg.at_node['topographic__steepest_slope']
    s18 = sa_factor * (a18 ** -0.5)
    s28 = sa_factor * (a28 ** -0.5)
    assert_equal(np.round(s[18], 3), np.round(s18, 3))
    assert_equal(np.round(s[28], 3), np.round(s28, 3))

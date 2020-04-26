#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Thu Jul 27 14:23:25 2017

@author: gtucker
"""

import numpy as np
from numpy.testing import assert_equal

from landlab import RasterModelGrid
from landlab.components import ErosionDeposition, FlowAccumulator


def test_erodep_slope_area_small_vs():
    """Test steady state run with Vs << 1."""

    # Set up a 5x5 grid with open boundaries and low initial elevations.
    rg = RasterModelGrid((5, 5))
    z = rg.add_zeros("topographic__elevation", at="node")
    z[:] = 0.01 * rg.x_of_node

    # Create a D8 flow handler
    fa = FlowAccumulator(rg, flow_director="FlowDirectorD8")

    # Parameter values for test 1
    K = 0.001
    vs = 0.0001
    U = 0.001
    dt = 10.0

    # Create the ErosionDeposition component...
    ed = ErosionDeposition(rg, K=K, v_s=vs, m_sp=0.5, n_sp=1.0, solver="adaptive")

    # ... and run it to steady state.
    for i in range(1000):
        fa.run_one_step()
        ed.run_one_step(dt=dt)
        z[rg.core_nodes] += U * dt

    # Test the results
    s = rg.at_node["topographic__steepest_slope"]
    sa_factor = (1.0 + vs) * U / K
    a11 = 2.0
    a12 = 1.0
    s = rg.at_node["topographic__steepest_slope"]
    s11 = sa_factor * (a11 ** -0.5)
    s12 = sa_factor * (a12 ** -0.5)
    assert_equal(np.round(s[11], 3), np.round(s11, 3))
    assert_equal(np.round(s[12], 3), np.round(s12, 3))


def test_erodep_slope_area_big_vs():
    """Test steady state run with Vs >> 1."""

    # Set up a 5x5 grid with open boundaries and low initial elevations.
    rg = RasterModelGrid((5, 5))
    z = rg.add_zeros("topographic__elevation", at="node")
    z[:] = 0.01 * rg.x_of_node

    # Create a D8 flow handler
    fa = FlowAccumulator(rg, flow_director="FlowDirectorD8")

    # Next test: big Vs
    K = 1.0
    vs = 1000.0
    U = 0.001
    dt = 10.0

    # Create the ErosionDeposition component...
    ed = ErosionDeposition(rg, K=K, v_s=vs, m_sp=0.5, n_sp=1.0, solver="adaptive")

    # ... and run it to steady state.
    for i in range(1000):
        fa.run_one_step()
        ed.run_one_step(dt=dt)
        z[rg.core_nodes] += U * dt

    # Test the results
    s = rg.at_node["topographic__steepest_slope"]
    sa_factor = (1.0 + vs) * U / K
    a11 = 2.0
    a12 = 1.0
    s11 = sa_factor * (a11 ** -0.5)
    s12 = sa_factor * (a12 ** -0.5)
    assert_equal(np.round(s[11], 2), np.round(s11, 2))
    assert_equal(np.round(s[12], 2), np.round(s12, 2))


def test_erodep_slope_area_with_vs_unity():
    """Test steady state run with Vs = 1."""

    # Set up a 5x5 grid with open boundaries and low initial elevations.
    rg = RasterModelGrid((5, 5))
    z = rg.add_zeros("topographic__elevation", at="node")
    z[:] = 0.01 * rg.x_of_node

    # Create a D8 flow handler
    fa = FlowAccumulator(rg, flow_director="FlowDirectorD8")

    # test: Vs = 1
    K = 0.002
    vs = 1.0
    U = 0.001
    dt = 10.0

    # Create the ErosionDeposition component...
    ed = ErosionDeposition(rg, K=K, v_s=vs, m_sp=0.5, n_sp=1.0, solver="adaptive")

    # ... and run it to steady state.
    for i in range(1000):
        fa.run_one_step()
        ed.run_one_step(dt=dt)
        z[rg.core_nodes] += U * dt

    # Test the results
    s = rg.at_node["topographic__steepest_slope"]
    sa_factor = (1.0 + vs) * U / K
    a11 = 2.0
    a12 = 1.0
    s11 = sa_factor * (a11 ** -0.5)
    s12 = sa_factor * (a12 ** -0.5)
    assert_equal(np.round(s[11], 2), np.round(s11, 2))
    assert_equal(np.round(s[12], 2), np.round(s12, 2))


def test_erodep_slope_area_shear_stress_scaling():
    """Test steady state run with m_sp = 0.33, n_sp=0.67, Vs = 1."""

    # Set up a 5x5 grid with open boundaries and low initial elevations.
    rg = RasterModelGrid((5, 5))
    rg.set_closed_boundaries_at_grid_edges(True, True, True, False)
    z = rg.add_zeros("topographic__elevation", at="node")
    z[:] = 0.01 * rg.x_of_node

    # Create a D8 flow handler
    fa = FlowAccumulator(rg, flow_director="FlowDirectorD8")

    # test: Vs = 1
    K = 0.002
    vs = 1.0
    U = 0.001
    dt = 10.0
    m_sp = 0.33
    n_sp = 0.67
    # Create the ErosionDeposition component...
    ed = ErosionDeposition(rg, K=K, v_s=vs, m_sp=m_sp, n_sp=n_sp, solver="adaptive")

    # ... and run it to steady state.
    for i in range(1500):
        fa.run_one_step()
        ed.run_one_step(dt=dt)
        z[rg.core_nodes] += U * dt

    # Test the results
    s = rg.at_node["topographic__steepest_slope"]
    sa_factor = ((1.0 + vs) * U / K) ** (1.0 / n_sp)
    a6 = rg.at_node["drainage_area"][6]
    a8 = rg.at_node["drainage_area"][8]
    s6 = sa_factor * (a6 ** -(m_sp / n_sp))
    s8 = sa_factor * (a8 ** -(m_sp / n_sp))
    assert_equal(np.round(s[6], 2), np.round(s6, 2))
    assert_equal(np.round(s[8], 2), np.round(s8, 2))


def test_erodep_slope_area_with_threshold():
    """Test steady state run with Vs = 1 and wc = 0.00001."""

    # Set up a 5x5 grid with open boundaries and low initial elevations.
    rg = RasterModelGrid((5, 5))
    z = rg.add_zeros("topographic__elevation", at="node")
    z[:] = 0.01 * rg.x_of_node

    # Create a D8 flow handler
    fa = FlowAccumulator(rg, flow_director="FlowDirectorD8")

    # test: Vs = 1
    K = 0.002
    vs = 1.0
    U = 0.001
    dt = 10.0
    wc = 0.0001

    # Create the ErosionDeposition component...
    ed = ErosionDeposition(
        rg, K=K, v_s=vs, m_sp=0.5, n_sp=1.0, sp_crit=wc, solver="adaptive"
    )

    # ... and run it to steady state.
    for i in range(1000):
        fa.run_one_step()
        ed.run_one_step(dt=dt)
        z[rg.core_nodes] += U * dt

    # Test the results
    s = rg.at_node["topographic__steepest_slope"]
    sa_factor = ((1.0 + vs) * U + wc) / K  # approximate sol'n
    a11 = 2.0
    a12 = 1.0
    s11 = sa_factor * (a11 ** -0.5)
    s12 = sa_factor * (a12 ** -0.5)
    assert_equal(np.round(s[11], 2), np.round(s11, 2))
    assert_equal(np.round(s[12], 2), np.round(s12, 2))

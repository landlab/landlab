#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Thu Jul 27 14:23:25 2017

@author: gtucker
"""

import numpy as np
import pytest
from numpy import testing

from landlab import HexModelGrid, RasterModelGrid
from landlab.components import ErosionDeposition, FlowAccumulator


def test_route_to_multiple_error_raised():
    mg = RasterModelGrid((10, 10))
    z = mg.add_zeros("topographic__elevation", at="node")
    z += mg.x_of_node + mg.y_of_node
    fa = FlowAccumulator(mg, flow_director="MFD")
    fa.run_one_step()

    with pytest.raises(NotImplementedError):
        ErosionDeposition(mg, K=0.01, v_s=0.001, m_sp=0.5, n_sp=1.0, sp_crit=0)


def test_phi_error_raised():
    mg = RasterModelGrid((10, 10))
    z = mg.add_zeros("topographic__elevation", at="node")
    z += mg.x_of_node + mg.y_of_node
    fa = FlowAccumulator(mg)
    fa.run_one_step()

    with pytest.raises(ValueError):
        ErosionDeposition(mg, phi=0)


def test_extra_kwd_error_raised():
    mg = RasterModelGrid((10, 10))
    z = mg.add_zeros("topographic__elevation", at="node")
    z += mg.x_of_node + mg.y_of_node
    fa = FlowAccumulator(mg)
    fa.run_one_step()

    with pytest.raises(ValueError):
        ErosionDeposition(mg, spam=0)


def test_bad_solver_name():
    """
    Test that any solver name besides 'basic' and 'adaptive' raises an error.
    """

    # set up a 5x5 grid with one open outlet node and low initial elevations.
    nr = 5
    nc = 5
    mg = RasterModelGrid((nr, nc), xy_spacing=10.0)

    mg.add_zeros("topographic__elevation", at="node")

    mg["node"]["topographic__elevation"] += (
        mg.node_y / 10000 + mg.node_x / 10000 + np.random.rand(len(mg.node_y)) / 10000
    )
    mg.set_closed_boundaries_at_grid_edges(
        bottom_is_closed=True,
        left_is_closed=True,
        right_is_closed=True,
        top_is_closed=True,
    )
    mg.set_watershed_boundary_condition_outlet_id(
        0, mg["node"]["topographic__elevation"], -9999.0
    )

    # Create a D8 flow handler
    FlowAccumulator(mg, flow_director="D8")

    # try to instantiate ErodionDeposition using a wrong solver name
    with pytest.raises(ValueError):
        ErosionDeposition(
            mg,
            K=0.01,
            v_s=0.001,
            m_sp=0.5,
            n_sp=1.0,
            sp_crit=0,
            F_f=0.0,
            solver="something_else",
        )


def test_steady_state_with_basic_solver_option():
    """
    Test that model matches the transport-limited analytical solution
    for slope/area relationship at steady state: S=((U * v_s) / (K * A^m)
    + U / (K * A^m))^(1/n).

    Also test that model matches the analytical solution for steady-state
    sediment flux: Qs = U * A * (1 - phi).
    """

    # set up a 5x5 grid with one open outlet node and low initial elevations.
    nr = 5
    nc = 5
    mg = RasterModelGrid((nr, nc), xy_spacing=10.0)

    z = mg.add_zeros("topographic__elevation", at="node")

    mg["node"]["topographic__elevation"] += (
        mg.node_y / 100000 + mg.node_x / 100000 + np.random.rand(len(mg.node_y)) / 10000
    )
    mg.set_closed_boundaries_at_grid_edges(
        bottom_is_closed=True,
        left_is_closed=True,
        right_is_closed=True,
        top_is_closed=True,
    )
    mg.set_watershed_boundary_condition_outlet_id(
        0, mg["node"]["topographic__elevation"], -9999.0
    )

    # Create a D8 flow handler
    fa = FlowAccumulator(
        mg, flow_director="D8", depression_finder="DepressionFinderAndRouter"
    )

    # Parameter values for detachment-limited test
    K = 0.01
    U = 0.0001
    dt = 1.0
    F_f = 0.0  # all sediment is considered coarse bedload
    m_sp = 0.5
    n_sp = 1.0
    v_s = 0.5

    # Instantiate the ErosionDeposition component...
    ed = ErosionDeposition(
        mg, K=K, F_f=F_f, v_s=v_s, m_sp=m_sp, n_sp=n_sp, sp_crit=0, solver="basic",
    )

    # ... and run it to steady state (5000x1-year timesteps).
    for i in range(5000):
        fa.run_one_step()
        ed.run_one_step(dt=dt)
        z[mg.core_nodes] += U * dt  # m

    # compare numerical and analytical slope solutions
    num_slope = mg.at_node["topographic__steepest_slope"][mg.core_nodes]
    analytical_slope = np.power(
        ((U * v_s) / (K * np.power(mg.at_node["drainage_area"][mg.core_nodes], m_sp)))
        + ((U) / (K * np.power(mg.at_node["drainage_area"][mg.core_nodes], m_sp))),
        1.0 / n_sp,
    )

    # test for match with analytical slope-area relationship
    testing.assert_array_almost_equal(
        num_slope,
        analytical_slope,
        decimal=8,
        err_msg="E/D slope-area test failed",
        verbose=True,
    )

    # compare numerical and analytical sediment flux solutions
    num_sedflux = mg.at_node["sediment__flux"][mg.core_nodes]
    analytical_sedflux = U * mg.at_node["drainage_area"][mg.core_nodes]

    # test for match with anakytical sediment flux
    testing.assert_array_almost_equal(
        num_sedflux,
        analytical_sedflux,
        decimal=8,
        err_msg="E/D sediment flux test failed",
        verbose=True,
    )


def test_can_run_with_hex():
    """Test that model can run with hex model grid."""

    # Set up a 5x5 grid with open boundaries and low initial elevations.
    mg = HexModelGrid((7, 7))
    z = mg.add_zeros("topographic__elevation", at="node")
    z[:] = 0.01 * mg.x_of_node

    # Create a D8 flow handler
    fa = FlowAccumulator(mg, flow_director="FlowDirectorSteepest")

    # Parameter values for test 1
    K = 0.001
    vs = 0.0001
    U = 0.001
    dt = 10.0

    # Create the ErosionDeposition component...
    ed = ErosionDeposition(mg, K=K, v_s=vs, m_sp=0.5, n_sp=1.0, solver="adaptive")

    # ... and run it to steady state.
    for i in range(2000):
        fa.run_one_step()
        ed.run_one_step(dt=dt)
        z[mg.core_nodes] += U * dt

    # Test the results
    s = mg.at_node["topographic__steepest_slope"]
    sa_factor = (1.0 + vs) * U / K
    a18 = mg.at_node["drainage_area"][18]
    a28 = mg.at_node["drainage_area"][28]
    s = mg.at_node["topographic__steepest_slope"]
    s18 = sa_factor * (a18 ** -0.5)
    s28 = sa_factor * (a28 ** -0.5)
    testing.assert_equal(np.round(s[18], 3), np.round(s18, 3))
    testing.assert_equal(np.round(s[28], 3), np.round(s28, 3))

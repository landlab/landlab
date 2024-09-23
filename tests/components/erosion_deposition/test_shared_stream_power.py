#!/usr/bin/env python2
"""
Created on Thu Jul 27 14:23:25 2017

@author: gtucker
"""

import numpy as np
import pytest
from numpy import testing

from landlab import HexModelGrid
from landlab import RasterModelGrid
from landlab.components import FlowAccumulator
from landlab.components import SharedStreamPower


def test_route_to_multiple_error_raised():
    grid = RasterModelGrid((10, 10))
    grid.at_node["topographic__elevation"] = grid.x_of_node + grid.y_of_node
    fa = FlowAccumulator(grid, flow_director="MFD")
    fa.run_one_step()

    with pytest.raises(NotImplementedError):
        SharedStreamPower(grid)


def test_bad_solver_name():
    """Test that any solver name besides 'basic' and 'adaptive' raises an error."""

    grid = RasterModelGrid((5, 5), xy_spacing=10.0)
    with pytest.raises(ValueError):
        SharedStreamPower(grid, solver="something_else")


def test_steady_state_with_basic_solver_option():
    """
    Test that model matches the transport-limited analytical solution
    for slope/area relationship at steady state: S=((U * v_s) / (K * A^m)
    + U / (K * A^m))^(1/n).

    Also test that model matches the analytical solution for steady-state
    sediment flux: Qs = U * A * (1 - phi).
    """

    # set up a 5x5 grid with one open outlet node and low initial elevations.
    mg = RasterModelGrid((5, 5), xy_spacing=10.0)

    mg.at_node["topographic__elevation"] = (
        mg.node_y / 100000 + mg.node_x / 100000 + np.random.rand(len(mg.node_y)) / 10000
    )
    mg.set_closed_boundaries_at_grid_edges(
        bottom_is_closed=True,
        left_is_closed=True,
        right_is_closed=True,
        top_is_closed=True,
    )
    mg.set_watershed_boundary_condition_outlet_id(
        0, mg.at_node["topographic__elevation"], -9999.0
    )

    # Create a D8 flow handler
    fa = FlowAccumulator(
        mg, flow_director="D8", depression_finder="DepressionFinderAndRouter"
    )

    # Instantiate the SharedStreamPower component...
    ed = SharedStreamPower(
        mg,
        k_bedrock=0.01,
        k_transport=0.01 / 0.05,
        F_f=0.0,  # all sediment is considered coarse bedload
        # v_s=v_s,
        m_sp=0.5,
        n_sp=1.0,
        sp_crit=0.0,
        solver="basic",
    )

    # ... and run it to steady state (5000x1-year timesteps).
    uplift = 0.0001
    dt = 1.0
    for _ in range(5000):
        fa.run_one_step()
        ed.run_one_step(dt=dt)
        mg.at_node["topographic__elevation"][mg.core_nodes] += uplift * dt

    analytical_slope = np.power(
        (
            (uplift * ed.v_s)
            / (
                ed.k_bedrock
                * np.power(mg.at_node["drainage_area"][mg.core_nodes], ed.m_sp)
            )
        )
        + (
            uplift
            / (
                ed.k_bedrock
                * np.power(mg.at_node["drainage_area"][mg.core_nodes], ed.m_sp)
            )
        ),
        1.0 / ed.n_sp,
    )

    # compare numerical and analytical slope solutions
    testing.assert_array_almost_equal(
        mg.at_node["topographic__steepest_slope"][mg.core_nodes],
        analytical_slope,
        decimal=8,
        err_msg="E/D slope-area test failed",
        verbose=True,
    )

    actual = mg.at_node["sediment__flux"][mg.core_nodes]
    expected = uplift * mg.at_node["drainage_area"][mg.core_nodes]

    # test for match with anakytical sediment flux
    testing.assert_array_almost_equal(
        actual,
        expected,
        decimal=8,
        err_msg="E/D sediment flux test failed",
        verbose=True,
    )


def test_can_run_with_hex():
    """Test that model can run with hex model grid."""

    # Set up a 5x5 grid with open boundaries and low initial elevations.
    mg = HexModelGrid((7, 7))
    mg.at_node["topographic__elevation"] = 0.01 * mg.x_of_node

    # Create a D8 flow handler
    fa = FlowAccumulator(mg, flow_director="FlowDirectorSteepest")

    # Create the SharedStreamPower component...
    ed = SharedStreamPower(
        mg,
        k_bedrock=0.001,
        k_transport=0.001 / 0.0001,
        m_sp=0.5,
        n_sp=1.0,
        solver="adaptive",
    )

    # ... and run it to steady state.
    uplift = 0.001
    dt = 10.0
    for _ in range(2000):
        fa.run_one_step()
        ed.run_one_step(dt=dt)
        mg.at_node["topographic__elevation"][mg.core_nodes] += uplift * dt

    sa_factor = (1.0 + ed.v_s) * uplift / ed.k_bedrock

    testing.assert_almost_equal(
        mg.at_node["topographic__steepest_slope"][mg.core_nodes],
        sa_factor * mg.at_node["drainage_area"][mg.core_nodes] ** -0.5,
    )

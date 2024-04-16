#! /usr/bin/env python
"""
Unit tests for landlab.components.gravel_bedrock_eroder.gravel_bedrock_eroder
"""
import numpy as np
from numpy.testing import assert_almost_equal

from landlab import HexModelGrid
from landlab.components import FlowAccumulator
from landlab.components import GravelBedrockEroder


def test_transport_rate():
    grid = HexModelGrid((4, 2), spacing=1000.0)
    grid.status_at_node[grid.perimeter_nodes] = grid.BC_NODE_IS_CLOSED
    grid.status_at_node[0] = grid.BC_NODE_IS_FIXED_VALUE

    elev = grid.add_zeros("topographic__elevation", at="node")
    elev[:] = (0.01 * grid.y_of_node) / np.cos(np.radians(30.0))
    sed = grid.add_zeros("soil__depth", at="node")
    sed[:] = 100.0

    fa = FlowAccumulator(grid)
    fa.run_one_step()
    gbe = GravelBedrockEroder(grid, intermittency_factor=0.02, depth_decay_scale=0.5)
    rock = grid.at_node["bedrock__elevation"]
    qs_out = grid.at_node["bedload_sediment__volume_outflux"]

    gbe.run_one_step(1.0e-6)  # using dt=0 prevents change to sed, rock, or elev
    assert_almost_equal(qs_out[grid.core_nodes], [9.88854526, 3.29618175, 3.29618175])

    elev[:] = (0.01 * grid.y_of_node) / np.cos(np.radians(30.0))
    sed[:] = 0.5
    rock[:] = elev - sed

    gbe.run_one_step(1.0e-6)
    assert_almost_equal(qs_out[grid.core_nodes], [6.25075275, 2.08358425, 2.08358425])

    elev[:] = (0.01 * grid.y_of_node) / np.cos(np.radians(30.0))
    sed[:] = 0.0
    rock[:] = elev

    gbe.run_one_step(1.0e-6)
    assert_almost_equal(qs_out[grid.core_nodes], [0.0, 0.0, 0.0])


def test_sediment_abrasion_rate():
    grid = HexModelGrid((4, 2), spacing=1000.0)
    grid.status_at_node[grid.perimeter_nodes] = grid.BC_NODE_IS_CLOSED
    grid.status_at_node[0] = grid.BC_NODE_IS_FIXED_VALUE

    elev = grid.add_zeros("topographic__elevation", at="node")
    elev[:] = (0.01 * grid.y_of_node) / np.cos(np.radians(30.0))
    sed = grid.add_zeros("soil__depth", at="node")
    sed[:] = 100.0

    fa = FlowAccumulator(grid)
    fa.run_one_step()
    gbe = GravelBedrockEroder(grid, abrasion_coefficient=1.0e-4)
    gbe.run_one_step(1.0)

    assert_almost_equal(
        gbe._abrasion[grid.core_nodes],
        [4.7576285545378313e-07, 9.515257103302159e-08, 9.515257103302159e-08],
    )


def test_rock_abrasion_rate():
    grid = HexModelGrid((4, 2), spacing=1000.0)
    grid.status_at_node[grid.perimeter_nodes] = grid.BC_NODE_IS_CLOSED
    grid.status_at_node[0] = grid.BC_NODE_IS_FIXED_VALUE

    elev = grid.add_zeros("topographic__elevation", at="node")
    elev[:] = (0.01 * grid.y_of_node) / np.cos(np.radians(30.0))
    sed = grid.add_zeros("soil__depth", at="node")
    sed[:] = 1.0

    fa = FlowAccumulator(grid)
    fa.run_one_step()
    gbe = GravelBedrockEroder(grid, abrasion_coefficient=1.0e-4)
    gbe.run_one_step(1.0)

    assert_almost_equal(
        gbe._sediment_outflux[grid.core_nodes], [3.12537638, 1.04179213, 1.04179213]
    )
    assert_almost_equal(
        gbe._rock_abrasion_rate[grid.core_nodes],
        [1.10635873e-07, 2.21271745e-08, 2.21271745e-08],
    )


def test_rock_plucking_rate():
    grid = HexModelGrid((4, 2), spacing=1000.0)
    grid.status_at_node[grid.perimeter_nodes] = grid.BC_NODE_IS_CLOSED
    grid.status_at_node[0] = grid.BC_NODE_IS_FIXED_VALUE

    elev = grid.add_zeros("topographic__elevation", at="node")
    elev[:] = (0.01 * grid.y_of_node) / np.cos(np.radians(30.0))
    sed = grid.add_zeros("soil__depth", at="node")
    sed[:] = 1.0

    fa = FlowAccumulator(grid)
    fa.run_one_step()
    gbe = GravelBedrockEroder(grid, plucking_coefficient=1.0e-4)
    gbe.run_one_step(1.0)

    assert_almost_equal(
        gbe._pluck_rate[grid.core_nodes],
        [5.12263532e-06, 1.70754511e-06, 1.70754511e-06],
    )


def test_steady_unlimited_sediment():
    grid = HexModelGrid((4, 2), spacing=1000.0)
    grid.status_at_node[grid.perimeter_nodes] = grid.BC_NODE_IS_CLOSED
    grid.status_at_node[0] = grid.BC_NODE_IS_FIXED_VALUE

    elev = grid.add_zeros("topographic__elevation", at="node")
    elev[:] = (0.13 * grid.y_of_node) / np.cos(np.radians(30.0))
    sed = grid.add_zeros("soil__depth", at="node")
    sed[:] = 10000.0
    rock = grid.add_zeros("bedrock__elevation", at="node")
    rock[:] = elev - sed

    fa = FlowAccumulator(grid)
    fa.run_one_step()
    gbe = GravelBedrockEroder(grid, abrasion_coefficient=0.0005)

    dt = 4.0e4
    uplift_rate = 0.0001
    nsteps = 500
    for _ in range(nsteps):
        elev[grid.core_nodes] += uplift_rate * dt
        sed[grid.core_nodes] += uplift_rate * dt
        gbe.run_one_step(dt)

    assert_almost_equal(
        grid.at_node["bedload_sediment__volume_outflux"][grid.core_nodes],
        [99.073, 45.033, 45.033],
        decimal=2,
    )
    assert_almost_equal(
        grid.at_node["topographic__steepest_slope"][grid.core_nodes],
        [0.130579, 0.170346, 0.170346],
        decimal=5,
    )


def test_steady_general():
    grid = HexModelGrid((3, 2), spacing=1000.0)
    grid.status_at_node[grid.perimeter_nodes] = grid.BC_NODE_IS_CLOSED
    grid.status_at_node[0] = grid.BC_NODE_IS_FIXED_VALUE

    elev = grid.add_zeros("topographic__elevation", at="node")
    elev[:] = (0.2 * grid.y_of_node) / np.cos(np.radians(30.0))
    sed = grid.add_zeros("soil__depth", at="node")
    sed[:] = 1.0
    rock = grid.add_zeros("bedrock__elevation", at="node")
    rock[:] = elev - sed

    fa = FlowAccumulator(grid)
    fa.run_one_step()
    gbe = GravelBedrockEroder(
        grid, abrasion_coefficient=0.0005, coarse_fraction_from_plucking=0.5
    )

    dt = 7500.0
    uplift_rate = 0.0001
    nsteps = 3000
    for _ in range(nsteps):
        elev[grid.core_nodes] += uplift_rate * dt
        rock[grid.core_nodes] += uplift_rate * dt
        gbe.run_one_step(dt)

    assert_almost_equal(
        grid.at_node["bedrock__exposure_fraction"][grid.core_nodes], 0.5062, decimal=4
    )
    assert_almost_equal(
        grid.at_node["topographic__steepest_slope"][grid.core_nodes], 0.2387, decimal=4
    )
    assert_almost_equal(
        grid.at_node["bedload_sediment__volume_outflux"][grid.core_nodes],
        32.972,
        decimal=3,
    )

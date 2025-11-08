#! /usr/bin/env python
"""
Unit tests for landlab.components.gravel_bedrock_eroder.gravel_bedrock_eroder
"""
import numpy as np
from numpy.testing import assert_almost_equal

from landlab import HexModelGrid
from landlab import RasterModelGrid
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
    gbe = GravelBedrockEroder(grid, abrasion_coefficients=[1.0e-4])
    gbe.run_one_step(1.0)

    assert_almost_equal(
        gbe._sed_abr_rates[0, grid.core_nodes],
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
    gbe = GravelBedrockEroder(grid, abrasion_coefficients=[1.0e-4])
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
    """
    Test the steady-state solution for the case of unlimited sediment.

    Notes
    -----
    The test case uses a HexModelGrid with 3 core nodes, two of which drain
    to the 3rd, which then drains to a single open boundary node.

    At steady state, the incoming sediment from uplift and (for the downstream
    node) upstream delivery must match the sum of outflux and loss to abrasion.
    With U as uplift rate, A as the surface area of a cell, Q discharge, S slope,
    Kqi as the product of transport coefficient and efficiency factor, p as the
    porosity, Qs as sediment flux, b as abrasion coefficient, and dx as
    link length,

    Qs_out + b (Qs_out + Qs_in) dx / 2 = U A (1 - p)

    Qs_out (1 + b dx / 2) = U A - b dx Qs_in / 2

    Qs_out = (U A (1 - p) - b dx Qs_in / 2) / (1 + b dx / 2)

    Kqi Q S^(7/6) = (U A (1 - p) - b dx Qs_in / 2) / (1 + b dx / 2)

    S = ((U A (1 - p) - b dx Qs_in / 2) / ((1 + b dx / 2) Kqi Q))^(6 / 7)

    For each of the 2 upstream nodes, Qs_in = 0, so

    Qs_out = U A (1 - p) / (1 + b dx / 2) ~ 0.0001 m/y 866,025 m2 / (1 + 0.0005 1000 / 2)
    ~ 45 m3/y

    S ~ 0.0001 m/y 866,025 m2 / ((1 + 0.0005 1000 / 2) 0.00041 866,025 m3/y
    ~0.170346

    (Calculation not shown for downstream nodes, but follows the same math)
    """
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
    gbe = GravelBedrockEroder(grid, abrasion_coefficients=[0.0005])

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
        decimal=0,
    )
    assert_almost_equal(
        grid.at_node["topographic__steepest_slope"][grid.core_nodes],
        [0.130579, 0.170346, 0.170346],
        decimal=3,
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
        grid, abrasion_coefficients=[0.0005], coarse_fractions_from_plucking=[0.5]
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


def test_calculations_in_sequence():
    """
    Test the results of calculations performed by run_one_step(), in sequence.

    run_one_step() calls update_rates(), then _update_rock_sed_and_elev()

    update_rates() updates:
        - slope gradient: here should be 0.001
        - rock exposure fraction: here should be ~0.5
        - tranport rate (outflux) at the core node: here should be
          0.05 x 0.01 x 1e6 m3/y x (0.001 m/m)^(7/6) x 0.5 ~ 0.079 m3/y
        - sediment influx: at node 1 (open boundary) should equal above
        - bedrock plucking rate: should be
          0.0001 x 0.01 x 1e6 x 0.001^(7/6) x 0.5 / 1e3 ~ 1.58113883e-7 m/y
        - sediment abrasion rate: should be
          0.5 x (0.079 + 0) m3/y x 0.0001 1/m x 1000 m / 1e6 m2 ~ 3.952847e-9 m/y
        - bedrock abrasion rate: should be 0.5 times above ~ 1.9764235e-09 m/y
        - sediment rate of change: should be
          (1 / (1 - 0.5)) x ((0 - 0.079) / 1e6) + 1.58e-7 - 3.95e-9)
          = 2 x (7.9e-8 + 1.58e-7 - 3.95e-9)
          ~ +1.501e-7
        - rock lowering rate: should be ~1.58e-7 + 2.0e-9 ~ 1.60e-7

    _update_rock_sed_and_elev() updates:
        - sediment thickness:
          starting with 0.69314718 m, adding 10,000 x 1.50208e-7 ~0.69464926 m
        - bedrock elevation:
          starting with 0.3068528 m, losing ~1.6e-7 x 10,000 ~0.3052528
    """
    grid = RasterModelGrid((3, 3), 1000.0)
    grid.status_at_node[grid.perimeter_nodes] = grid.BC_NODE_IS_CLOSED
    grid.status_at_node[1] = grid.BC_NODE_IS_FIXED_VALUE

    elev = grid.add_zeros("topographic__elevation", at="node")
    elev[4] = 1.0
    sed = grid.add_zeros("soil__depth", at="node")
    sed[:] = -np.log(0.5)
    rock = grid.add_zeros("bedrock__elevation", at="node")
    rock[:] = elev - sed

    fa = FlowAccumulator(grid)
    fa.run_one_step()
    gbe = GravelBedrockEroder(
        grid,
        transport_coefficient=0.05,
        sediment_porosity=0.5,
        depth_decay_scale=1.0,
        plucking_coefficient=1.0e-4,
        number_of_sediment_classes=1,
        abrasion_coefficients=[0.0001],
        coarse_fractions_from_plucking=[1.0],
    )
    gbe.run_one_step(
        0.001
    )  # tiny time step to calculate rates but mostly preserve elev and sed thickness

    assert_almost_equal(grid.at_node["bedrock__elevation"][4], 0.30685281944, decimal=6)
    assert_almost_equal(grid.at_node["soil__depth"][4], 0.69314718, decimal=6)
    assert_almost_equal(grid.at_node["topographic__steepest_slope"][4], 0.001)
    assert_almost_equal(grid.at_node["bedrock__exposure_fraction"][4], 0.5)
    assert_almost_equal(
        grid.at_node["bedload_sediment__volume_outflux"][4], 0.0790569415
    )
    assert_almost_equal(
        grid.at_node["bedload_sediment__volume_influx"][1], 0.0790569415
    )
    assert_almost_equal(
        grid.at_node["bedrock__plucking_rate"][4], 1.58113883e-7, decimal=15
    )
    print("passed to here, trying sar")
    assert_almost_equal(
        gbe._sed_abr_rates[0, 4],
        3.952847e-9,
        decimal=15,
    )
    assert_almost_equal(
        grid.at_node["bedrock__abrasion_rate"][4], 1.9764235e-09, decimal=15
    )
    assert_almost_equal(gbe._dHdt[4], 1.50e-7, decimal=9)
    assert_almost_equal(gbe._rock_lowering_rate[4], 1.60e-7, decimal=9)

    gbe.run_one_step(
        10000.0
    )  # big time step to get measurable change in elev and sed thickness

    assert_almost_equal(sed[4], 0.694649262)
    assert_almost_equal(rock[4], 0.305252, decimal=6)
    assert_almost_equal(elev[4], 0.9999, decimal=5)

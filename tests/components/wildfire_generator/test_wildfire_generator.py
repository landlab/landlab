#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jun  5 16:24:02 2026

@author: matheuswanderleydealmeida

Unit tests for wildfire generation using the WildfireGenerator component
"""

import numpy as np
import pytest
from numpy import testing

from landlab import RasterModelGrid
from landlab import FieldError
from landlab.components import WildfireGenerator


def test_missing_fuel_availability():
    """
    WildfireGenerator should raise a FieldError when the fuel_availability
    field is not present on the grid.
    """
    mg = RasterModelGrid((5, 5), xy_spacing=1)
    mg.add_zeros("topographic__elevation", at="node")
    mg.add_zeros("drainage_area", at="node")

    with pytest.raises(FieldError):
        WildfireGenerator(mg)


def test_fire_frequency():
    """
    WildfireGenerator should produce a fire_log of length 17 after the
    first run_one_step call with a fixed seed.
    """
    np.random.seed(5000)
    dt = 1
    mg = RasterModelGrid((5, 5), xy_spacing=1.0)
    mg.add_zeros("topographic__elevation", at="node")
    mg.add_zeros("drainage_area", at="node")
    fuel = mg.add_zeros("fuel_availability", at="node")
    fuel[:] = np.random.rand(mg.number_of_nodes) * 0.75

    wg = WildfireGenerator(mg)
    wg.run_one_step(dt)

    testing.assert_equal(len(wg.fire_log), 17)


def test_fire_locations():
    """
    WildfireGenerator should ignite fires at the expected node coordinates
    given a fixed random seed.
    """
    np.random.seed(5000)
    dt = 1
    mg = RasterModelGrid((5, 5), xy_spacing=1.0)
    mg.add_zeros("topographic__elevation", at="node")
    mg.add_zeros("drainage_area", at="node")
    fuel = mg.add_zeros("fuel_availability", at="node")
    fuel[:] = np.random.rand(mg.number_of_nodes) * 0.75
    wg = WildfireGenerator(mg)

    wg.run_one_step(dt)

    expected_locations = [
        np.array([3., 1.]), np.array([1., 1.]), np.array([2., 1.]),
        np.array([3., 2.]), np.array([1., 3.]), np.array([1., 2.]),
        np.array([1., 1.]), np.array([2., 2.]), np.array([1., 2.]),
        np.array([2., 3.]), np.array([2., 2.]), np.array([2., 3.]),
        np.array([2., 3.]), np.array([3., 2.]), np.array([2., 2.]),
        np.array([3., 3.]), np.array([1., 1.]),
    ]

    actual_locations = wg.ignition_nodes

    assert len(actual_locations) == len(expected_locations)
    for actual, expected in zip(actual_locations, expected_locations):
        testing.assert_array_equal(actual, expected)


def test_regrow_vegetation_changing_dt():
    """
    _regrow_vegetation should scale with dt according to the recovery
    timescales from _get_regrowth_params.
    """
    mg = RasterModelGrid((5, 5), xy_spacing=1.0)
    mg.add_zeros("topographic__elevation", at="node")
    mg.add_zeros("drainage_area", at="node")
    mg.add_zeros("fuel_availability", at="node")

    wg = WildfireGenerator(mg)
    regrowth_time, regrowth_exponent = wg._get_regrowth_params()

    for dt in [0.5, 1, 5, 20]:
        wg._vegetation[:] = 0.0
        wg._regrow_vegetation(dt)
        expected = 1 - np.exp(-regrowth_exponent * dt / regrowth_time)
        testing.assert_allclose(wg._vegetation, expected)


def test_missing_drainage_area_raises():
    """
    WildfireGenerator should raise a FieldError when the drainage_area
    field is not present on the grid.
    """
    mg = RasterModelGrid((5, 5), xy_spacing=1)
    mg.add_zeros("topographic__elevation", at="node")
    mg.add_ones("fuel_availability", at="node")

    with pytest.raises(FieldError):
        WildfireGenerator(mg)


def test_vegetation_field_linked_to_grid():
    """
    The _vegetation array inside WildfireGenerator should reference the same
    memory as the grid's fuel_availability field.
    """
    mg = RasterModelGrid((9, 9), xy_spacing=1)
    z = mg.add_zeros("topographic__elevation", at="node")
    z[:] = 100
    mg.add_ones("fuel_availability", at="node")
    mg.add_zeros("drainage_area", at="node")
    mg.at_node["fuel_availability"][:] = 0.6
    wf = WildfireGenerator(mg)

    testing.assert_array_equal(wf._vegetation, mg.at_node["fuel_availability"])


def test_rivers_mask_based_on_minimum_river_threshold():
    """
    Nodes with drainage_area > minimum_river_threshold should be marked
    as rivers in the _rivers property.
    """
    mg = RasterModelGrid((9, 9), xy_spacing=1)
    z = mg.add_zeros("topographic__elevation", at="node")
    z[:] = 100
    mg.add_ones("fuel_availability", at="node")
    mg.add_zeros("drainage_area", at="node")
    mg.at_node["drainage_area"][20] = 3e5   # above default threshold of 2e5
    mg.at_node["drainage_area"][21] = 1e5   # below

    wf = WildfireGenerator(mg)

    assert wf._rivers[20] == True
    assert wf._rivers[21] == False


@pytest.mark.parametrize("aridity, expected_t, expected_e", [
    (0.05, 5.4,  0.42),   # Trees + shrubs
    (0.2,  8.0,  0.29),   # Dense shrubs
    (0.4,  12.4, 0.18),   # Dense shrubs high end
    (0.6,  15.0, 0.15),   # Sparse shrubs
    (0.9,  18.5, 0.12),   # Sparse shrubs high end
    (1.0,  18.5, 0.12),   # Boundary value aridity == 1
])
def test_get_regrowth_params(aridity, expected_t, expected_e):
    """
    _get_regrowth_params should return the correct (t, e) pair for each
    aridity bin, including the boundary value aridity == 1.
    """
    mg = RasterModelGrid((9, 9), xy_spacing=1)
    z = mg.add_zeros("topographic__elevation", at="node")
    z[:] = 100
    mg.add_ones("fuel_availability", at="node")
    mg.add_zeros("drainage_area", at="node")
    wf = WildfireGenerator(mg, aridity=aridity)
    t, e = wf._get_regrowth_params()

    assert t == expected_t
    assert e == expected_e


@pytest.mark.parametrize("aridity, exp_fw, exp_aw", [
    (0.05, 0.5, 0.5),
    (0.2,  0.5, 0.5),
    (0.6,  0.8, 0.2),
    (0.9,  0.9, 0.1),
    (1.0,  0.9, 0.1),
])
def test_get_aridity_weights(aridity, exp_fw, exp_aw):
    """
    _get_aridity_weights should return the correct (fuel_weight, aridity_weight)
    pair for each aridity bin.
    """
    mg = RasterModelGrid((9, 9), xy_spacing=1)
    z = mg.add_zeros("topographic__elevation", at="node")
    z[:] = 100
    mg.add_ones("fuel_availability", at="node")
    mg.add_zeros("drainage_area", at="node")
    wf = WildfireGenerator(mg, aridity=aridity)
    fw, aw = wf._get_aridity_weights()

    assert fw == exp_fw
    assert aw == exp_aw


def test_calc_severity_range():
    """
    _calc_severity should always return a value in [0, 1].
    """
    mg = RasterModelGrid((9, 9), xy_spacing=1)
    z = mg.add_zeros("topographic__elevation", at="node")
    z[:] = 100
    mg.add_ones("fuel_availability", at="node")
    mg.add_zeros("drainage_area", at="node")
    wf = WildfireGenerator(mg, aridity=0.5)

    for _ in range(200):
        s = wf._calc_severity()
        assert 0.0 <= s <= 1.0, f"Severity {s} is out of [0, 1]"


def test_calc_severity_rounded():
    """
    _calc_severity should return a value rounded to 2 decimal places.
    """
    mg = RasterModelGrid((9, 9), xy_spacing=1)
    z = mg.add_zeros("topographic__elevation", at="node")
    z[:] = 100
    mg.add_ones("fuel_availability", at="node")
    mg.add_zeros("drainage_area", at="node")
    wf = WildfireGenerator(mg)

    for _ in range(50):
        s = wf._calc_severity()
        assert round(s, 2) == s, f"Severity {s} is not rounded to 2 decimals"


def test_calc_severity_increases_with_aridity():
    """
    On average, higher aridity should produce higher severity (power-law
    relationship). Test using a large sample to reduce noise.
    """
    mg_low = RasterModelGrid((9, 9), xy_spacing=1)
    z = mg_low.add_zeros("topographic__elevation", at="node")
    z[:] = 100
    mg_low.add_ones("fuel_availability", at="node")
    mg_low.add_zeros("drainage_area", at="node")

    mg_high = RasterModelGrid((9, 9), xy_spacing=1)
    z = mg_high.add_zeros("topographic__elevation", at="node")
    z[:] = 100
    mg_high.add_ones("fuel_availability", at="node")
    mg_high.add_zeros("drainage_area", at="node")

    wf_low = WildfireGenerator(mg_low, aridity=0.1)
    wf_high = WildfireGenerator(mg_high, aridity=0.9)

    np.random.seed(0)
    severities_low = [wf_low._calc_severity() for _ in range(500)]
    severities_high = [wf_high._calc_severity() for _ in range(500)]

    assert np.mean(severities_high) > np.mean(severities_low), (
        "Higher aridity should produce higher mean severity"
    )


def test_regrow_vegetation_increases_values():
    """
    _regrow_vegetation should increase vegetation values at all nodes
    (when not already at max).
    """
    mg = RasterModelGrid((9, 9), xy_spacing=1)
    z = mg.add_zeros("topographic__elevation", at="node")
    z[:] = 100
    mg.add_ones("fuel_availability", at="node")
    mg.add_zeros("drainage_area", at="node")
    mg.at_node["fuel_availability"][:] = 0.5
    wf = WildfireGenerator(mg, aridity=0.5)

    veg_before = wf._vegetation.copy()
    wf._regrow_vegetation(dt=1)

    assert np.all(wf._vegetation >= veg_before), (
        "Vegetation should not decrease after regrowth"
    )


def test_regrow_vegetation_capped_at_max():
    """
    _regrow_vegetation should never push vegetation above max_vegetation.
    """
    mg = RasterModelGrid((9, 9), xy_spacing=1)
    z = mg.add_zeros("topographic__elevation", at="node")
    z[:] = 100
    mg.add_ones("fuel_availability", at="node")
    mg.add_zeros("drainage_area", at="node")
    mg.at_node["fuel_availability"][:] = 1.0   # already at max
    wf = WildfireGenerator(mg, aridity=0.5)

    wf._regrow_vegetation(dt=1)

    assert np.all(wf._vegetation <= wf._max_vegetation), (
        "Vegetation exceeded max_vegetation after regrowth"
    )


def test_regrow_vegetation_recovers_over_time():
    """
    After starting at zero vegetation, repeated regrowth steps should
    eventually bring vegetation back close to max_vegetation.
    """
    mg = RasterModelGrid((9, 9), xy_spacing=1)
    z = mg.add_zeros("topographic__elevation", at="node")
    z[:] = 100
    mg.add_zeros("fuel_availability", at="node")
    mg.add_zeros("drainage_area", at="node")
    wf = WildfireGenerator(mg, aridity=0.05)   # fast regrowth bin

    for _ in range(500):
        wf._regrow_vegetation(dt=1)

    assert np.all(wf._vegetation > 0.9), (
        "Vegetation should recover close to max_vegetation after many regrowth steps"
    )


def test_get_neighbors_returns_valid_nodes():
    """
    _get_neighbors should return only valid (non-negative) active node indices
    for an interior node.
    """
    mg = RasterModelGrid((9, 9), xy_spacing=1)
    z = mg.add_zeros("topographic__elevation", at="node")
    z[:] = 100
    mg.add_ones("fuel_availability", at="node")
    mg.add_zeros("drainage_area", at="node")
    wf = WildfireGenerator(mg)

    interior_node = 40   # centre of a 9x9 grid
    neighbors = wf._get_neighbors(interior_node)

    assert len(neighbors) > 0, "Interior node should have neighbors"
    for n in neighbors:
        assert n >= 0, f"Neighbor index {n} is invalid (boundary placeholder)"


def test_get_neighbors_interior_node_count():
    """
    A fully interior node on a 9x9 grid should have exactly 4 D4 neighbors.
    """
    mg = RasterModelGrid((9, 9), xy_spacing=1)
    z = mg.add_zeros("topographic__elevation", at="node")
    z[:] = 100
    mg.add_ones("fuel_availability", at="node")
    mg.add_zeros("drainage_area", at="node")
    wf = WildfireGenerator(mg)

    neighbors = wf._get_neighbors(40)   # centre node, fully surrounded
    active_neighbors = [n for n in neighbors if n >= 0]

    assert len(active_neighbors) == 4


def test_fire_reduces_vegetation():
    """
    Nodes recorded in fire_log "burned_nodes" should have lower vegetation
    than the initial value of 1.
    """
    np.random.seed(7)
    mg = RasterModelGrid((9, 9), xy_spacing=1)
    z = mg.add_zeros("topographic__elevation", at="node")
    z[:] = 100
    mg.add_ones("fuel_availability", at="node")
    mg.add_zeros("drainage_area", at="node")
    wf = WildfireGenerator(mg, potential_fires=200, aridity=0.9)
    wf._current_time = 1

    wf._fire(dt=1)

    if len(wf.fire_log) > 0:
        burned_nodes = wf.fire_log[0]["burned_nodes"]
        assert np.all(wf._vegetation[burned_nodes] < 1.0), (
            "Vegetation at burned nodes should be reduced below 1"
        )


def test_fire_does_not_cross_rivers():
    """
    Fire should not burn nodes flagged as rivers.
    """
    np.random.seed(0)
    mg = RasterModelGrid((9, 9), xy_spacing=1)
    z = mg.add_zeros("topographic__elevation", at="node")
    z[:] = 100
    mg.add_ones("fuel_availability", at="node")
    mg.add_zeros("drainage_area", at="node")

    # Mark a column of nodes as rivers (high drainage area)
    river_nodes = [4, 13, 22, 31, 40, 49, 58, 67, 76]
    mg.at_node["drainage_area"][river_nodes] = 3e5

    wf = WildfireGenerator(mg, potential_fires=300, aridity=0.9)

    for t in range(20):
        wf._current_time = t
        wf._fire(dt=1)

    if len(wf.fire_log) > 0:
        all_burned = set(
            n for record in wf.fire_log for n in record["burned_nodes"]
        )
        for rn in river_nodes:
            assert rn not in all_burned, (
                f"River node {rn} should not be burned"
            )


# =============================================================================
# run_one_step tests
# =============================================================================

def test_run_one_step_updates_vegetation():
    """
    run_one_step should modify the vegetation field (fire + regrowth combined)
    while keeping it within [0, max_vegetation].
    """
    np.random.seed(42)
    mg = RasterModelGrid((9, 9), xy_spacing=1)
    z = mg.add_zeros("topographic__elevation", at="node")
    z[:] = 100
    mg.add_ones("fuel_availability", at="node")
    mg.add_zeros("drainage_area", at="node")
    mg.at_node["fuel_availability"][:] = 0.5
    wf = WildfireGenerator(mg, potential_fires=100, aridity=0.6)

    wf.run_one_step(dt=1)

    assert np.all(wf._vegetation >= 0), "Vegetation dropped below 0"
    assert np.all(wf._vegetation <= wf._max_vegetation), "Vegetation exceeded max"


def test_run_one_step_appends_fire_log():
    """
    Repeated calls to run_one_step should accumulate entries in fire_log.
    """
    np.random.seed(5)
    mg = RasterModelGrid((9, 9), xy_spacing=1)
    z = mg.add_zeros("topographic__elevation", at="node")
    z[:] = 100
    mg.add_ones("fuel_availability", at="node")
    mg.add_zeros("drainage_area", at="node")
    wf = WildfireGenerator(mg, potential_fires=200, aridity=0.9)

    for _ in range(10):
        wf.run_one_step(dt=1)

    assert len(wf.fire_log) >= 1, (
        "fire_log should accumulate entries over multiple steps"
    )
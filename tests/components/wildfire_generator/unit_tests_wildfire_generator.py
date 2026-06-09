#!/usr/bin/env python3
"""
Created on Fri Jun  5 16:24:02 2026

@author: matheuswanderleydealmeida

Unit tests for wildfire generation using the WildfireGenerator component
"""

import numpy as np
import pandas as pd
import pytest
from numpy import testing

from landlab import RasterModelGrid
from landlab.components import WildfireGenerator


def make_grid(n_rows=9, n_columns=9, spacing=1, plateau_elev=100.0):
    """Helper to create a standard RasterModelGrid for testing."""
    mg = RasterModelGrid((n_rows, n_columns), xy_spacing=spacing)
    z = mg.add_zeros("topographic__elevation", at="node")
    z[:] = plateau_elev
    mg.add_ones("fuel_availability", at="node")
    mg.add_zeros("drainage_area", at="node")
    return mg


def test_missing_fuel_availability_raises():
    """
    WildfireGenerator should raise a KeyError when the fuel_availability
    field is not present on the grid.
    """
    mg = RasterModelGrid((5, 5), xy_spacing=1)
    mg.add_zeros("topographic__elevation", at="node")
    mg.add_zeros("drainage_area", at="node")
    # mg.add_zeros("fuel_availability", at="node")

    with pytest.raises(KeyError):
        WildfireGenerator(mg)


def test_missing_drainage_area_raises():
    """
    WildfireGenerator should raise a KeyError when the drainage_area
    field is not present on the grid.
    """
    mg = RasterModelGrid((5, 5), xy_spacing=1)
    mg.add_zeros("topographic__elevation", at="node")
    mg.add_ones("fuel_availability", at="node")
    # mg.add_zeros("drainage_area", at="node")

    with pytest.raises(KeyError):
        WildfireGenerator(mg)


def test_instantiation_default_params():
    """
    WildfireGenerator should instantiate correctly with a valid grid and
    default parameters, and expose the expected attributes.
    """
    mg = make_grid()
    wf = WildfireGenerator(mg)

    assert wf.potential_fires == 100
    assert wf.dt == 1
    assert wf.dx == 1
    assert wf.alpha == 0.3
    assert wf.aridity == 0.5
    assert wf.sev_exponent == 0.64
    assert wf.last_fire_time is None
    assert isinstance(wf.fire_log, pd.DataFrame)
    assert len(wf.fire_log) == 0


def test_instantiation_custom_params():
    """
    WildfireGenerator should store custom parameters correctly.
    """
    mg = make_grid()
    wf = WildfireGenerator(mg, potential_fires=50, aridity=0.8, alpha=0.3)

    assert wf.potential_fires == 50
    assert wf.aridity == 0.8
    assert wf.alpha == 0.3


def test_vegetation_field_linked_to_grid():
    """
    The vegetation array inside WildfireGenerator should reference the same
    memory as the grid field, so changes to one are reflected in the other.
    """
    mg = make_grid()
    mg.at_node["fuel_availability"][:] = 0.6
    wf = WildfireGenerator(mg)

    testing.assert_array_equal(wf.vegetation, mg.at_node["fuel_availability"])


def test_rivers_mask_based_on_riv_min():
    """
    Nodes with drainage_area > riv_min should be marked as rivers.
    """
    mg = make_grid()
    mg.at_node["drainage_area"][20] = 3e5  # above default riv_min of 2e5
    mg.at_node["drainage_area"][21] = 1e5  # below

    wf = WildfireGenerator(mg)

    assert wf.rivers[20] is True
    assert wf.rivers[21] is False


@pytest.mark.parametrize(
    "aridity, expected_t, expected_e",
    [
        (0.05, 5.4, 0.42),  # Trees + shrubs
        (0.2, 8.0, 0.29),  # Dense shrubs
        (0.4, 12.4, 0.18),  # Dense shrubs high end
        (0.6, 15.0, 0.15),  # Sparse shrubs
        (0.9, 18.5, 0.12),  # Sparse shrubs high end
        (1.0, 18.5, 0.12),  # Boundary value aridity == 1
    ],
)
def test_get_regrowth_params(aridity, expected_t, expected_e):
    """
    _get_regrowth_params should return the correct (t, e) pair for each
    aridity bin, including the boundary value aridity == 1.
    """
    mg = make_grid()
    wf = WildfireGenerator(mg, aridity=aridity)
    t, e = wf._get_regrowth_params()

    assert t == expected_t
    assert e == expected_e


@pytest.mark.parametrize(
    "aridity, exp_fw, exp_aw",
    [
        (0.05, 0.5, 0.5),
        (0.2, 0.5, 0.5),
        (0.6, 0.8, 0.2),
        (0.9, 0.9, 0.1),
        (1.0, 0.9, 0.1),
    ],
)
def test_get_aridity_weights(aridity, exp_fw, exp_aw):
    """
    get_aridity_weights should return the correct (fuel_weight, aridity_weight)
    pair for each aridity bin.
    """
    mg = make_grid()
    wf = WildfireGenerator(mg, aridity=aridity)
    fw, aw = wf.get_aridity_weights()

    assert fw == exp_fw
    assert aw == exp_aw


def test_calc_severity_range():
    """
    calc_severity should always return a value in [0, 1].
    """
    mg = make_grid()
    wf = WildfireGenerator(mg, aridity=0.5)

    for _ in range(200):
        s = wf.calc_severity()
        assert 0.0 <= s <= 1.0, f"Severity {s} is out of [0, 1]"


def test_calc_severity_rounded():
    """
    calc_severity should return a value rounded to 2 decimal places.
    """
    mg = make_grid()
    wf = WildfireGenerator(mg)

    for _ in range(50):
        s = wf.calc_severity()
        assert round(s, 2) == s, f"Severity {s} is not rounded to 2 decimals"


def test_calc_severity_increases_with_aridity():
    """
    On average, higher aridity should produce higher severity (power-law
    relationship). Test using a large sample to reduce noise.
    """
    mg_low = make_grid()
    mg_high = make_grid()

    wf_low = WildfireGenerator(mg_low, aridity=0.1)
    wf_high = WildfireGenerator(mg_high, aridity=0.9)

    np.random.seed(0)
    severities_low = [wf_low.calc_severity() for _ in range(500)]
    severities_high = [wf_high.calc_severity() for _ in range(500)]

    assert np.mean(severities_high) > np.mean(
        severities_low
    ), "Higher aridity should produce higher mean severity"


def test_regrow_vegetation_increases_values():
    """
    regrow_vegetation should increase vegetation values at all nodes.
    """
    mg = make_grid()
    mg.at_node["fuel_availability"][:] = 0.5
    wf = WildfireGenerator(mg, aridity=0.5)

    veg_before = wf.vegetation.copy()
    wf.regrow_vegetation(t=1)

    assert np.all(
        wf.vegetation >= veg_before
    ), "Vegetation should not decrease after regrowth"


def test_regrow_vegetation_capped_at_max():
    """
    regrow_vegetation should never push vegetation above max_vegetation (1.0).
    """
    mg = make_grid()
    mg.at_node["fuel_availability"][:] = 1.0  # already at max
    wf = WildfireGenerator(mg, aridity=0.5)

    wf.regrow_vegetation(t=1)

    assert np.all(
        wf.vegetation <= wf.max_vegetation
    ), "Vegetation exceeded max_vegetation after regrowth"


def test_regrow_vegetation_recovers_over_time():
    """
    After a complete burn, repeated regrowth steps should eventually bring
    vegetation back close to max_vegetation.
    """
    mg = make_grid()
    mg.at_node["fuel_availability"][:] = 0.0
    wf = WildfireGenerator(mg, aridity=0.05)  # fast regrowth bin

    for t in range(500):
        wf.regrow_vegetation(t)

    assert np.all(
        wf.vegetation > 0.9
    ), "Vegetation should recover close to 1 after many regrowth steps"


def test_get_neighbors_returns_valid_nodes():
    """
    get_neighbors should return only valid (non-negative) active node indices.
    """
    mg = make_grid()
    wf = WildfireGenerator(mg)

    interior_node = 40  # centre of a 9x9 grid
    neighbors = wf.get_neighbors(interior_node)

    assert len(neighbors) > 0, "Interior node should have neighbors"
    for n in neighbors:
        assert n >= 0, f"Neighbor index {n} is invalid (boundary placeholder)"


def test_get_neighbors_interior_node_count():
    """
    A fully interior node on a 9x9 grid should have exactly 4 D4 neighbors.
    """
    mg = make_grid()
    wf = WildfireGenerator(mg)

    neighbors = wf.get_neighbors(40)  # centre node, fully surrounded
    active_neighbors = [n for n in neighbors if n >= 0]

    assert len(active_neighbors) == 4


def test_fire_returns_list():
    """
    fire() should return a list (possibly empty).
    """
    mg = make_grid()
    wf = WildfireGenerator(mg, potential_fires=10, aridity=0.5)

    result = wf.fire(t=1)
    assert isinstance(result, (list, set)), "fire() should return a list or set"


def test_fire_log_populated_after_fire():
    """
    After a successful fire event, fire_log should contain at least one row
    with the expected columns.
    """
    np.random.seed(42)
    mg = make_grid()
    mg.at_node["fuel_availability"][:] = 1.0
    wf = WildfireGenerator(mg, potential_fires=200, aridity=0.9)

    wf.fire(t=5)

    expected_cols = {
        "year",
        "fire_size (km2)",
        "center_X",
        "center_Y",
        "severity_factor",
        "aridity",
        "changed_nodes",
        "pct of vegetation removed",
    }
    assert expected_cols.issubset(
        set(wf.fire_log.columns)
    ), "fire_log is missing expected columns"
    assert len(wf.fire_log) >= 1, "fire_log should have at least one entry"
    assert (wf.fire_log["year"] == 5).all(), "fire_log year should match t"


def test_fire_reduces_vegetation():
    """
    Nodes recorded in fire_log changed_nodes should have lower vegetation
    than the initial value of 1.
    """
    np.random.seed(7)
    mg = make_grid()
    mg.at_node["fuel_availability"][:] = 1.0
    wf = WildfireGenerator(mg, potential_fires=200, aridity=0.9)

    wf.fire(t=1)

    if len(wf.fire_log) > 0:
        burned_nodes = wf.fire_log["changed_nodes"].iloc[0]
        assert np.all(
            wf.vegetation[burned_nodes] < 1.0
        ), "Vegetation at burned nodes should be reduced below 1"


def test_fire_does_not_cross_rivers():
    """
    Fire should not burn nodes flagged as rivers.
    """
    np.random.seed(0)
    mg = make_grid()
    mg.at_node["fuel_availability"][:] = 1.0

    # Mark a column of nodes as rivers (high drainage area)
    river_nodes = [4, 13, 22, 31, 40, 49, 58, 67, 76]
    mg.at_node["drainage_area"][river_nodes] = 3e5

    wf = WildfireGenerator(mg, potential_fires=300, aridity=0.9)

    for t in range(20):
        wf.fire(t)

    if len(wf.fire_log) > 0:
        all_burned = {n for nodes in wf.fire_log["changed_nodes"] for n in nodes}
        for rn in river_nodes:
            assert rn not in all_burned, f"River node {rn} should not be burned"


def test_fire_initialises_last_fire_time():
    """
    After the first call to fire(), last_fire_time should be initialised as
    an array of the correct length.
    """
    np.random.seed(1)
    mg = make_grid()
    mg.at_node["fuel_availability"][:] = 1.0
    wf = WildfireGenerator(mg, potential_fires=50)

    assert wf.last_fire_time is None
    wf.fire(t=3)
    assert wf.last_fire_time is not None
    assert len(wf.last_fire_time) == mg.number_of_nodes


# =============================================================================
# run_one_step tests
# =============================================================================


def test_run_one_step_updates_vegetation():
    """
    run_one_step should modify the vegetation field (fire + regrowth combined).
    """
    np.random.seed(42)
    mg = make_grid()
    mg.at_node["fuel_availability"][:] = 0.5
    wf = WildfireGenerator(mg, potential_fires=100, aridity=0.6)

    wf.run_one_step(t=1)

    # After fire+regrowth, vegetation values may go up or down, but they must
    # remain in [0, max_vegetation].
    assert np.all(wf.vegetation >= 0), "Vegetation dropped below 0"
    assert np.all(wf.vegetation <= wf.max_vegetation), "Vegetation exceeded max"


def test_run_one_step_appends_fire_log():
    """
    Repeated calls to run_one_step should accumulate entries in fire_log.
    """
    np.random.seed(5)
    mg = make_grid()
    mg.at_node["fuel_availability"][:] = 1.0
    wf = WildfireGenerator(mg, potential_fires=200, aridity=0.9)

    for t in range(10):
        wf.run_one_step(t)

    assert (
        len(wf.fire_log) >= 1
    ), "fire_log should accumulate entries over multiple steps"

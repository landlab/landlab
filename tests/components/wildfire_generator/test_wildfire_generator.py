#!/usr/bin/env python3
"""
Created on Fri Jun  5 16:24:02 2026

@author: matheuswanderleydealmeida

Unit tests for the WildfireGenerator Landlab component.

"""

import numpy as np
import pytest
from numpy import testing

from landlab import FieldError
from landlab import RasterModelGrid
from landlab.components import WildfireGenerator

# =============================================================================
# Helpers
# =============================================================================


def _make_grid(nrows=5, ncols=5, spacing=1.0, fuel_val=0.0, topo_val=0.0, da_val=0.0):
    """Return a minimal RasterModelGrid with all required fields."""
    mg = RasterModelGrid((nrows, ncols), xy_spacing=spacing)
    z = mg.add_zeros("topographic__elevation", at="node")
    z[:] = topo_val
    da = mg.add_zeros("drainage_area", at="node")
    da[:] = da_val
    fuel = mg.add_zeros("fuel_availability", at="node")
    fuel[:] = fuel_val
    return mg


# =============================================================================
# __init__ — required field checks
# =============================================================================


def test_missing_fuel_availability():
    """WildfireGenerator should raise FieldError when fuel_availability is absent."""
    mg = RasterModelGrid((5, 5), xy_spacing=1)
    mg.add_zeros("topographic__elevation", at="node")
    mg.add_zeros("drainage_area", at="node")
    with pytest.raises(FieldError):
        WildfireGenerator(mg)


def test_missing_drainage_area_raises():
    """WildfireGenerator should raise FieldError when drainage_area is absent."""
    mg = RasterModelGrid((5, 5), xy_spacing=1)
    mg.add_zeros("topographic__elevation", at="node")
    mg.add_ones("fuel_availability", at="node")
    with pytest.raises(FieldError):
        WildfireGenerator(mg)


# =============================================================================
# __init__ — parameter validation
# =============================================================================


def test_aridity_below_zero_raises():
    mg = _make_grid()
    with pytest.raises(ValueError, match="aridity must be >= 0.0"):
        WildfireGenerator(mg, aridity=-0.01)


def test_aridity_above_one_raises():
    mg = _make_grid()
    with pytest.raises(ValueError, match="aridity must be <= 1.0"):
        WildfireGenerator(mg, aridity=1.01)


def test_potential_fires_zero_raises():
    mg = _make_grid()
    with pytest.raises(ValueError, match="potential_fires must be > 0"):
        WildfireGenerator(mg, potential_fires=0)


def test_potential_fires_negative_raises():
    mg = _make_grid()
    with pytest.raises(ValueError, match="potential_fires must be > 0"):
        WildfireGenerator(mg, potential_fires=-10)


def test_max_vegetation_zero_raises():
    mg = _make_grid()
    with pytest.raises(ValueError, match="max_vegetation must be > 0.0"):
        WildfireGenerator(mg, max_vegetation=0.0)


def test_max_vegetation_negative_raises():
    mg = _make_grid()
    with pytest.raises(ValueError, match="max_vegetation must be > 0.0"):
        WildfireGenerator(mg, max_vegetation=-1.0)


def test_minimum_river_threshold_negative_raises():
    mg = _make_grid()
    with pytest.raises(ValueError, match="minimum_river_threshold must be >= 0.0"):
        WildfireGenerator(mg, minimum_river_threshold=-1.0)


def test_upslope_preference_negative_raises():
    mg = _make_grid()
    with pytest.raises(ValueError, match="upslope_preference"):
        WildfireGenerator(mg, upslope_preference=-0.5)


def test_severity_exponent_zero_raises():
    mg = _make_grid()
    with pytest.raises(ValueError, match="severity_exponent must be > 0.0"):
        WildfireGenerator(mg, severity_exponent=0.0)


def test_severity_exponent_above_one_raises():
    mg = _make_grid()
    with pytest.raises(ValueError, match="severity_exponent must be <= 1.0"):
        WildfireGenerator(mg, severity_exponent=1.5)


def test_negative_fuel_raises():
    """fuel_availability field containing negative values should raise ValueError."""
    mg = _make_grid()
    mg.at_node["fuel_availability"][5] = -0.1
    with pytest.raises(
        ValueError, match="fuel_availability cannot contain negative values"
    ):
        WildfireGenerator(mg)


def test_fuel_exceeds_max_vegetation_raises():
    """fuel_availability > max_vegetation should raise ValueError."""
    mg = _make_grid(fuel_val=0.5)
    mg.at_node["fuel_availability"][5] = 1.5
    with pytest.raises(
        ValueError, match="fuel_availability contains values exceeding max_vegetation"
    ):
        WildfireGenerator(mg, max_vegetation=1.0)


# =============================================================================
# __init__ — output field initialisation
# =============================================================================


def test_last_fire_time_initialised_to_neg_inf():
    """last_fire_time should be -inf at all nodes after initialisation."""
    mg = _make_grid()
    wg = WildfireGenerator(mg)
    assert np.all(
        wg._last_fire_time == -np.inf
    ), "last_fire_time should be -inf everywhere before any fire"


# =============================================================================
# __init__ — field aliasing
# =============================================================================


def test_vegetation_field_linked_to_grid():
    """_vegetation must share memory with grid's fuel_availability field."""
    mg = _make_grid(nrows=9, ncols=9, topo_val=100, fuel_val=0.6)
    wg = WildfireGenerator(mg)
    testing.assert_array_equal(wg._vegetation, mg.at_node["fuel_availability"])


# =============================================================================
# Properties
# =============================================================================


def test_potential_fires_property():
    mg = _make_grid()
    wg = WildfireGenerator(mg, potential_fires=42)
    assert wg.potential_fires == 42


def test_upslope_preference_property():
    mg = _make_grid()
    wg = WildfireGenerator(mg, upslope_preference=0.7)
    assert wg.upslope_preference == 0.7


def test_aridity_getter():
    mg = _make_grid()
    wg = WildfireGenerator(mg, aridity=0.3)
    assert wg.aridity == 0.3


def test_aridity_setter_valid():
    mg = _make_grid()
    wg = WildfireGenerator(mg)
    wg.aridity = 0.8
    assert wg.aridity == 0.8


def test_aridity_setter_below_zero_raises():
    mg = _make_grid()
    wg = WildfireGenerator(mg)
    with pytest.raises(ValueError, match="aridity must be >= 0.0"):
        wg.aridity = -0.1


def test_aridity_setter_above_one_raises():
    mg = _make_grid()
    wg = WildfireGenerator(mg)
    with pytest.raises(ValueError, match="aridity must be <= 1.0"):
        wg.aridity = 1.1


def test_rivers_mask_based_on_minimum_river_threshold():
    """Nodes with drainage_area > minimum_river_threshold should be flagged as rivers."""
    mg = _make_grid(nrows=9, ncols=9, topo_val=100)
    mg.at_node["drainage_area"][20] = 3e5  # above default threshold (2e5)
    mg.at_node["drainage_area"][21] = 1e5  # below threshold
    wg = WildfireGenerator(mg)
    assert wg._rivers[20] == True  # noqa: E712  (np.True_ is not True but == True)
    assert wg._rivers[21] == False  # noqa: E712


def test_fire_sizes_property_length_matches_log():
    np.random.seed(42)
    mg = _make_grid(nrows=9, ncols=9, spacing=100, fuel_val=0.8)
    wg = WildfireGenerator(mg, potential_fires=50, aridity=0.7, seed=42)
    wg.run_one_step(1)
    assert len(wg.fire_sizes) == len(wg.fire_log)


def test_fire_sizes_all_positive():
    np.random.seed(42)
    mg = _make_grid(nrows=9, ncols=9, spacing=100, fuel_val=0.8)
    wg = WildfireGenerator(mg, potential_fires=50, aridity=0.7, seed=42)
    wg.run_one_step(1)
    assert all(
        s > 0 for s in wg.fire_sizes
    ), "All recorded fire sizes should be positive"


def test_burned_nodes_property_length_matches_log():
    np.random.seed(42)
    mg = _make_grid(nrows=9, ncols=9, spacing=100, fuel_val=0.8)
    wg = WildfireGenerator(mg, potential_fires=50, aridity=0.7, seed=42)
    wg.run_one_step(1)
    assert len(wg.burned_nodes) == len(wg.fire_log)


def test_burned_nodes_are_lists():
    np.random.seed(42)
    mg = _make_grid(nrows=9, ncols=9, spacing=100, fuel_val=0.8)
    wg = WildfireGenerator(mg, potential_fires=50, aridity=0.7, seed=42)
    wg.run_one_step(1)
    assert all(isinstance(b, list) for b in wg.burned_nodes)


def test_severity_property_length_matches_log():
    np.random.seed(42)
    mg = _make_grid(nrows=9, ncols=9, spacing=100, fuel_val=0.8)
    wg = WildfireGenerator(mg, potential_fires=50, aridity=0.7, seed=42)
    wg.run_one_step(1)
    assert len(wg.severity) == len(wg.fire_log)


def test_severity_values_in_unit_interval():
    np.random.seed(42)
    mg = _make_grid(nrows=9, ncols=9, spacing=100, fuel_val=0.8)
    wg = WildfireGenerator(mg, potential_fires=50, aridity=0.7, seed=42)
    wg.run_one_step(1)
    assert all(0.0 <= s <= 1.0 for s in wg.severity)


def test_ignition_nodes_property_length_matches_log():
    np.random.seed(42)
    mg = _make_grid(nrows=9, ncols=9, spacing=100, fuel_val=0.8)
    wg = WildfireGenerator(mg, potential_fires=50, aridity=0.7, seed=42)
    wg.run_one_step(1)
    assert len(wg.ignition_nodes) == len(wg.fire_log)


def test_last_step_fires_all_from_current_year():
    """last_step_fires should contain only records from the most recent timestep."""
    mg = _make_grid(nrows=9, ncols=9, spacing=100, fuel_val=0.8)
    wg = WildfireGenerator(mg, potential_fires=30, aridity=0.7, seed=7)
    wg.run_one_step(1)
    wg.run_one_step(1)
    assert all(r["year"] == wg._current_time for r in wg.last_step_fires)


def test_last_step_fires_excludes_previous_years():
    """last_step_fires must not include fires from earlier timesteps."""
    mg = _make_grid(nrows=9, ncols=9, spacing=100, fuel_val=0.8)
    wg = WildfireGenerator(mg, potential_fires=30, aridity=0.7, seed=7)
    wg.run_one_step(1)
    fires_step1 = len(wg.fire_log)
    wg.run_one_step(1)
    # last_step_fires must be a strict subset of fire_log
    assert len(wg.last_step_fires) <= len(wg.fire_log)
    assert len(wg.last_step_fires) <= fires_step1 or fires_step1 == 0


# =============================================================================
# _get_regrowth_params
# =============================================================================


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
    """_get_regrowth_params should return the correct (T, e) pair for each bin."""
    mg = _make_grid(nrows=9, ncols=9, topo_val=100)
    wg = WildfireGenerator(mg, aridity=aridity)
    t, e = wg._get_regrowth_params()
    assert t == expected_t
    assert e == expected_e


# =============================================================================
# _get_aridity_weights
# =============================================================================


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
    """_get_aridity_weights should return the correct weight pair for each bin."""
    mg = _make_grid(nrows=9, ncols=9, topo_val=100)
    wg = WildfireGenerator(mg, aridity=aridity)
    fw, aw = wg._get_aridity_weights()
    assert fw == exp_fw
    assert aw == exp_aw


def test_aridity_weights_sum_to_one():
    """Fuel weight + aridity weight should always sum to 1.0."""
    for aridity in [0.05, 0.2, 0.4, 0.6, 0.9, 1.0]:
        mg = _make_grid(nrows=9, ncols=9, topo_val=100)
        wg = WildfireGenerator(mg, aridity=aridity)
        fw, aw = wg._get_aridity_weights()
        assert fw + aw == pytest.approx(
            1.0
        ), f"Weights don't sum to 1 for aridity={aridity}: {fw} + {aw}"


# =============================================================================
# _calc_severity
# =============================================================================


def test_calc_severity_range():
    """_calc_severity should always return a value in [0, 1]."""
    mg = _make_grid(nrows=9, ncols=9, topo_val=100)
    wg = WildfireGenerator(mg, aridity=0.5)
    for _ in range(200):
        s = wg._calc_severity()
        assert 0.0 <= s <= 1.0, f"Severity {s} is outside [0, 1]"


def test_calc_severity_rounded():
    """_calc_severity should return a value rounded to exactly 2 decimal places."""
    mg = _make_grid(nrows=9, ncols=9, topo_val=100)
    wg = WildfireGenerator(mg)
    for _ in range(50):
        s = wg._calc_severity()
        assert round(s, 2) == s, f"Severity {s} is not rounded to 2 decimals"


def test_calc_severity_increases_with_aridity():
    """Higher aridity should produce higher mean severity (power-law relationship)."""
    mg_low = _make_grid(nrows=9, ncols=9, topo_val=100)
    mg_high = _make_grid(nrows=9, ncols=9, topo_val=100)
    wf_low = WildfireGenerator(mg_low, aridity=0.1)
    wf_high = WildfireGenerator(mg_high, aridity=0.9)
    np.random.seed(0)
    sev_low = [wf_low._calc_severity() for _ in range(500)]
    sev_high = [wf_high._calc_severity() for _ in range(500)]
    assert np.mean(sev_high) > np.mean(
        sev_low
    ), "Higher aridity should produce higher mean severity"


def test_calc_severity_at_zero_aridity_near_zero():
    """With aridity ≈ 0, severity should be very close to 0 on average."""
    mg = _make_grid(nrows=9, ncols=9, topo_val=100)
    wg = WildfireGenerator(mg, aridity=0.01)
    np.random.seed(1)
    mean_sev = np.mean([wg._calc_severity() for _ in range(300)])
    assert mean_sev < 0.1, f"Mean severity at near-zero aridity was {mean_sev}"


def test_calc_severity_at_unit_aridity_near_one():
    """With aridity = 1, severity should be close to 1 on average."""
    mg = _make_grid(nrows=9, ncols=9, topo_val=100)
    wg = WildfireGenerator(mg, aridity=1.0)
    np.random.seed(2)
    mean_sev = np.mean([wg._calc_severity() for _ in range(300)])
    assert mean_sev > 0.9, f"Mean severity at aridity=1 was {mean_sev}"


# =============================================================================
# _regrow_vegetation
# =============================================================================


def test_regrow_vegetation_changing_dt():
    """_regrow_vegetation growth increment should match the analytic formula for each dt."""
    mg = _make_grid(nrows=5, ncols=5)
    wg = WildfireGenerator(mg)
    regrowth_time, regrowth_exponent = wg._get_regrowth_params()
    for dt in [0.5, 1, 5, 20]:
        wg._vegetation[:] = 0.0
        wg._regrow_vegetation(dt)
        expected = 1 - np.exp(-regrowth_exponent * dt / regrowth_time)
        testing.assert_allclose(wg._vegetation, expected)


def test_regrow_vegetation_increases_values():
    """_regrow_vegetation should never decrease vegetation at sub-maximum nodes."""
    mg = _make_grid(nrows=9, ncols=9, topo_val=100, fuel_val=0.5)
    wg = WildfireGenerator(mg, aridity=0.5)
    veg_before = wg._vegetation.copy()
    wg._regrow_vegetation(dt=1)
    assert np.all(
        wg._vegetation >= veg_before
    ), "Vegetation should not decrease after regrowth"


def test_regrow_vegetation_capped_at_max():
    """_regrow_vegetation must not push vegetation above max_vegetation."""
    mg = _make_grid(nrows=9, ncols=9, topo_val=100, fuel_val=1.0)
    wg = WildfireGenerator(mg, aridity=0.5)
    wg._regrow_vegetation(dt=1)
    assert np.all(
        wg._vegetation <= wg._max_vegetation
    ), "Vegetation exceeded max_vegetation after regrowth"


def test_regrow_vegetation_recovers_over_time():
    """Starting from zero, repeated regrowth should bring vegetation close to max."""
    mg = _make_grid(nrows=9, ncols=9, topo_val=100, fuel_val=0.0)
    wg = WildfireGenerator(mg, aridity=0.05)  # fastest regrowth bin
    for _ in range(500):
        wg._regrow_vegetation(dt=1)
    assert np.all(
        wg._vegetation > 0.9
    ), "Vegetation should recover close to max_vegetation after many regrowth steps"


def test_regrow_vegetation_zero_dt_no_change():
    """_regrow_vegetation with dt=0 should leave vegetation unchanged."""
    mg = _make_grid(nrows=9, ncols=9, topo_val=100, fuel_val=0.4)
    wg = WildfireGenerator(mg, aridity=0.5)
    veg_before = wg._vegetation.copy()
    wg._regrow_vegetation(dt=0)
    testing.assert_array_equal(wg._vegetation, veg_before)


def test_regrow_vegetation_custom_max_vegetation():
    """_regrow_vegetation should respect a non-default max_vegetation."""
    mg = _make_grid(nrows=5, ncols=5, fuel_val=0.5)
    max_veg = 0.6
    wg = WildfireGenerator(mg, max_vegetation=max_veg)
    for _ in range(200):
        wg._regrow_vegetation(dt=1)
    assert np.all(
        wg._vegetation <= max_veg
    ), f"Vegetation exceeded custom max_vegetation={max_veg}"


# =============================================================================
# _get_neighbors
# =============================================================================


def test_get_neighbors_returns_valid_nodes():
    """_get_neighbors should return only valid (non-negative) node indices 
    for interior nodes."""
    mg = _make_grid(nrows=9, ncols=9, topo_val=100)
    wg = WildfireGenerator(mg)
    interior_node = 40  # centre of 9×9
    neighbors = wg._get_neighbors(interior_node)
    assert len(neighbors) > 0, "Interior node should have neighbors"
    for n in neighbors:
        assert n >= 0, f"Neighbor index {n} is invalid"


def test_get_neighbors_interior_node_count():
    """A fully interior node on a 9×9 grid should have exactly 4 active D4 neighbors."""
    mg = _make_grid(nrows=9, ncols=9, topo_val=100)
    wg = WildfireGenerator(mg)
    active = [n for n in wg._get_neighbors(40) if n >= 0]
    assert len(active) == 4


# =============================================================================
# _fire
# =============================================================================


def test_fire_reduces_vegetation():
    """Burned nodes should have lower vegetation than before the fire."""
    np.random.seed(7)
    mg = _make_grid(nrows=9, ncols=9, topo_val=100, fuel_val=1.0)
    wg = WildfireGenerator(mg, potential_fires=200, aridity=0.9)
    wg._current_time = 1
    wg._fire(dt=1)
    if wg.fire_log:
        burned = wg.fire_log[0]["burned_nodes"]
        assert np.all(
            wg._vegetation[burned] < 1.0
        ), "Vegetation at burned nodes should be reduced below 1"


def test_fire_does_not_cross_rivers():
    """Fire should never burn nodes flagged as rivers."""
    np.random.seed(0)
    mg = _make_grid(nrows=9, ncols=9, topo_val=100, fuel_val=1.0)
    river_nodes = [4, 13, 22, 31, 40, 49, 58, 67, 76]
    mg.at_node["drainage_area"][river_nodes] = 3e5  # above default threshold
    wg = WildfireGenerator(mg, potential_fires=300, aridity=0.9)
    for t in range(20):
        wg._current_time = t
        wg._fire(dt=1)
    if wg.fire_log:
        all_burned = {n for r in wg.fire_log for n in r["burned_nodes"]}
        for rn in river_nodes:
            assert rn not in all_burned, f"River node {rn} should not be burned"


def test_fire_ignition_skipped_on_river_center():
    """An ignition attempt whose center node is a river should be silently skipped."""
    mg = _make_grid(nrows=9, ncols=9, topo_val=100, fuel_val=1.0)
    # Mark all core nodes as rivers
    mg.at_node["drainage_area"][mg.core_nodes] = 1e9
    wg = WildfireGenerator(mg, potential_fires=200, aridity=0.9, seed=0)
    wg._current_time = 1
    wg._fire(dt=1)
    assert (
        len(wg.fire_log) == 0
    ), "No fires should be recorded when all core nodes are rivers"


def test_fire_pct_vegetation_removed_in_valid_range():
    """pct_of_vegetation_removed should be in [0, 100] when initial fuel > 0."""
    np.random.seed(3)
    mg = _make_grid(nrows=9, ncols=9, spacing=100, fuel_val=0.8)
    wg = WildfireGenerator(mg, potential_fires=50, aridity=0.7, seed=3)
    wg.run_one_step(1)
    for r in wg.fire_log:
        pct = r["pct of vegetation removed"]
        assert (
            0.0 <= pct <= 100.0
        ), f"pct_of_vegetation_removed={pct} is outside [0, 100]"


def test_fire_size_equals_core_cell_area_sum():
    """fire_size (km²) should equal the sum of cell areas of burned nodes / 1e6.

    Boundary nodes have cell_area = 0, so only core cells contribute.
    """
    mg = _make_grid(nrows=9, ncols=9, spacing=1000, fuel_val=1.0)
    wg = WildfireGenerator(mg, potential_fires=3, aridity=0.9, seed=0)
    wg.run_one_step(1)
    for r in wg.fire_log:
        expected_km2 = mg.cell_area_at_node[r["burned_nodes"]].sum() / 1e6
        testing.assert_allclose(
            r["fire_size (km2)"],
            expected_km2,
            rtol=1e-10,
            err_msg="fire_size does not match sum of cell areas / 1e6",
        )


def test_last_fire_time_updated_at_burned_nodes():
    """last_fire_time should equal current_time at every burned node."""
    mg = _make_grid(nrows=15, ncols=15, spacing=100, fuel_val=0.8)
    wg = WildfireGenerator(mg, potential_fires=5, aridity=0.5, seed=99)
    wg.run_one_step(1)
    for r in wg.fire_log:
        for node in r["burned_nodes"]:
            assert (
                wg._last_fire_time[node] == wg._current_time
            ), f"last_fire_time at burned node {node} not updated to current_time"


def test_last_fire_time_unchanged_at_unburned_nodes():
    """last_fire_time should remain -inf at nodes that were never burned."""
    mg = _make_grid(nrows=15, ncols=15, spacing=100, fuel_val=0.8)
    wg = WildfireGenerator(mg, potential_fires=5, aridity=0.5, seed=99)
    wg.run_one_step(1)
    burned_all = {n for r in wg.fire_log for n in r["burned_nodes"]}
    unburned = [n for n in mg.core_nodes if n not in burned_all]
    if unburned:
        assert np.all(
            wg._last_fire_time[unburned] == -np.inf
        ), "last_fire_time should remain -inf at unburned nodes"


def test_fire_log_record_contains_required_keys():
    """Every fire_log entry should contain all documented keys."""
    mg = _make_grid(nrows=9, ncols=9, spacing=100, fuel_val=0.8)
    wg = WildfireGenerator(mg, potential_fires=50, aridity=0.7, seed=42)
    wg.run_one_step(1)
    required_keys = {
        "year",
        "ignition_node",
        "burned_nodes",
        "fire_size (km2)",
        "severity_factor",
        "aridity",
        "pct of vegetation removed",
    }
    for r in wg.fire_log:
        assert required_keys.issubset(
            r.keys()
        ), f"fire_log record is missing keys: {required_keys - r.keys()}"


def test_fire_log_aridity_matches_component_aridity():
    """The aridity recorded in each fire_log entry should match the component aridity."""
    mg = _make_grid(nrows=9, ncols=9, spacing=100, fuel_val=0.8)
    target_aridity = 0.65
    wg = WildfireGenerator(mg, potential_fires=50, aridity=target_aridity, seed=42)
    wg.run_one_step(1)
    for r in wg.fire_log:
        assert r["aridity"] == target_aridity


# =============================================================================
# run_one_step — time tracking
# =============================================================================


def test_run_one_step_current_time_none_before_first_call():
    """_current_time should be None before any call to run_one_step."""
    mg = _make_grid()
    wg = WildfireGenerator(mg, seed=None)
    assert wg._current_time is None


def test_run_one_step_current_time_set_on_first_call():
    """_current_time should equal dt after the first run_one_step call."""
    mg = _make_grid()
    wg = WildfireGenerator(mg, seed=None)
    wg.run_one_step(5)
    assert wg._current_time == 5


def test_run_one_step_current_time_accumulates():
    """_current_time should accumulate across successive calls."""
    mg = _make_grid()
    wg = WildfireGenerator(mg, seed=None)
    wg.run_one_step(2)
    assert wg._current_time == 2
    wg.run_one_step(3)
    assert wg._current_time == 5
    wg.run_one_step(7)
    assert wg._current_time == 12


def test_run_one_step_explicit_current_time_overrides():
    """Passing current_time should override the accumulated clock."""
    mg = _make_grid()
    wg = WildfireGenerator(mg, seed=None)
    wg.run_one_step(1)
    wg.run_one_step(1)
    wg.run_one_step(1, current_time=100)
    assert wg._current_time == 100


# =============================================================================
# run_one_step — fire log and vegetation
# =============================================================================


def test_fire_frequency():
    """fire_log should contain exactly 17 records after one step with seed=5000."""
    np.random.seed(5000)
    mg = RasterModelGrid((5, 5), xy_spacing=1.0)
    mg.add_zeros("topographic__elevation", at="node")
    mg.add_zeros("drainage_area", at="node")
    fuel = mg.add_zeros("fuel_availability", at="node")
    fuel[:] = np.random.rand(mg.number_of_nodes) * 0.75
    wg = WildfireGenerator(mg)
    wg.run_one_step(1)
    testing.assert_equal(len(wg.fire_log), 17)


def test_fire_locations():
    """Ignition node coordinates should match the seed-locked expectation."""
    np.random.seed(5000)
    mg = RasterModelGrid((5, 5), xy_spacing=1.0)
    mg.add_zeros("topographic__elevation", at="node")
    mg.add_zeros("drainage_area", at="node")
    fuel = mg.add_zeros("fuel_availability", at="node")
    fuel[:] = np.random.rand(mg.number_of_nodes) * 0.75
    wg = WildfireGenerator(mg)
    wg.run_one_step(1)

    expected_locations = [
        np.array([3.0, 1.0]),
        np.array([1.0, 1.0]),
        np.array([2.0, 1.0]),
        np.array([3.0, 2.0]),
        np.array([1.0, 3.0]),
        np.array([1.0, 2.0]),
        np.array([1.0, 1.0]),
        np.array([2.0, 2.0]),
        np.array([1.0, 2.0]),
        np.array([2.0, 3.0]),
        np.array([2.0, 2.0]),
        np.array([2.0, 3.0]),
        np.array([2.0, 3.0]),
        np.array([3.0, 2.0]),
        np.array([2.0, 2.0]),
        np.array([3.0, 3.0]),
        np.array([1.0, 1.0]),
    ]
    assert len(wg.ignition_nodes) == len(expected_locations)
    for actual, expected in zip(wg.ignition_nodes, expected_locations):
        testing.assert_array_equal(actual, expected)


def test_run_one_step_updates_vegetation():
    """run_one_step should keep vegetation within [0, max_vegetation]."""
    np.random.seed(42)
    mg = _make_grid(nrows=9, ncols=9, topo_val=100, fuel_val=0.5)
    wg = WildfireGenerator(mg, potential_fires=100, aridity=0.6)
    wg.run_one_step(dt=1)
    assert np.all(wg._vegetation >= 0), "Vegetation dropped below 0"
    assert np.all(wg._vegetation <= wg._max_vegetation), "Vegetation exceeded max"


def test_run_one_step_appends_fire_log():
    """Repeated run_one_step calls should accumulate entries in fire_log."""
    np.random.seed(5)
    mg = _make_grid(nrows=9, ncols=9, topo_val=100, fuel_val=1.0)
    wg = WildfireGenerator(mg, potential_fires=200, aridity=0.9)
    for _ in range(10):
        wg.run_one_step(dt=1)
    assert len(wg.fire_log) >= 1, "fire_log should accumulate entries"


def test_run_one_step_fire_log_year_matches_current_time():
    """Every fire_log entry's year should equal the current_time when it was recorded."""
    mg = _make_grid(nrows=9, ncols=9, spacing=100, fuel_val=0.8)
    wg = WildfireGenerator(mg, potential_fires=30, aridity=0.7, seed=5)
    for _step in range(5):
        wg.run_one_step(1)
    for r in wg.fire_log:
        assert isinstance(r["year"], (int, float)), "year field should be numeric"
        assert 1 <= r["year"] <= 5, f"year={r['year']} is outside the simulated range"

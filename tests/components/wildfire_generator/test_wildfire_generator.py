#!/usr/bin/env python3
"""
Created on Fri Jun  5 16:24:02 2026

@author: matheuswanderleydealmeida

Unit tests for the WildfireGenerator Landlab component.

"""

import numpy as np
import pytest
from numpy import testing
from requireit import ValidationError

from landlab import FieldError
from landlab import RasterModelGrid
from landlab.components import WildfireGenerator
from landlab.components.wildfire_generator import FireEvent
from landlab.components.wildfire_generator import FireEventLog
from landlab.components.wildfire_generator.wildfire_generator import get_aridity_weights
from landlab.components.wildfire_generator.wildfire_generator import get_regrowth_params

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


def _run_with_log(mg, n_steps=1, dt=1, **kwargs):
    """Run a WildfireGenerator with a FireEventLog and return both."""
    log = FireEventLog()
    wg = WildfireGenerator(mg, log=log, **kwargs)
    for _ in range(n_steps):
        wg.run_one_step(dt)
    return wg, log


def _doctest_grid():
    """Grid used by the class docstring example."""
    np.random.seed(5000)
    mg = RasterModelGrid((5, 5), xy_spacing=1.0)
    mg.add_zeros("topographic__elevation", at="node")
    mg.add_zeros("drainage_area", at="node")
    fuel = mg.add_zeros("fuel_availability", at="node")
    fuel[:] = np.random.rand(mg.number_of_nodes) * 0.75
    return mg


class _FakeRng:
    """Deterministic stand-in for numpy.random.Generator."""

    def poisson(self, lam):
        return 1

    def choice(self, arr):
        return 0

    def random(self):
        return 0.0

    def normal(self, loc, scale):
        return loc


# =============================================================================
# __init__ -- required field checks
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
# __init__ -- parameter validation
# =============================================================================


def test_aridity_below_zero_raises():
    mg = _make_grid()
    with pytest.raises(ValidationError, match="aridity must be >= 0.0"):
        WildfireGenerator(mg, aridity=-0.01)


def test_aridity_above_one_raises():
    mg = _make_grid()
    with pytest.raises(ValidationError, match="aridity must be <= 1.0"):
        WildfireGenerator(mg, aridity=1.01)


def test_potential_fires_zero_raises():
    mg = _make_grid()
    with pytest.raises(ValidationError, match="potential_fires must be > 0"):
        WildfireGenerator(mg, potential_fires=0)


def test_potential_fires_negative_raises():
    mg = _make_grid()
    with pytest.raises(ValidationError, match="potential_fires must be > 0"):
        WildfireGenerator(mg, potential_fires=-10)


def test_max_vegetation_zero_raises():
    mg = _make_grid()
    with pytest.raises(ValidationError, match="max_vegetation must be > 0.0"):
        WildfireGenerator(mg, max_vegetation=0.0)


def test_max_vegetation_negative_raises():
    mg = _make_grid()
    with pytest.raises(ValidationError, match="max_vegetation must be > 0.0"):
        WildfireGenerator(mg, max_vegetation=-1.0)


def test_minimum_river_threshold_zero_raises():
    mg = _make_grid()
    with pytest.raises(ValidationError, match="minimum_river_threshold must be > 0.0"):
        WildfireGenerator(mg, minimum_river_threshold=0.0)


def test_minimum_river_threshold_negative_raises():
    mg = _make_grid()
    with pytest.raises(ValidationError, match="minimum_river_threshold must be > 0.0"):
        WildfireGenerator(mg, minimum_river_threshold=-1.0)


def test_minimum_river_threshold_inf_allowed():
    """np.inf disables firebreaks and should be accepted."""
    mg = _make_grid()
    WildfireGenerator(mg, minimum_river_threshold=np.inf)


def test_upslope_preference_negative_raises():
    mg = _make_grid()
    with pytest.raises(ValidationError, match="upslope_preference"):
        WildfireGenerator(mg, upslope_preference=-0.5)


def test_severity_exponent_negative_raises():
    mg = _make_grid()
    with pytest.raises(ValidationError, match="severity_exponent must be >= 0.0"):
        WildfireGenerator(mg, severity_exponent=-0.1)


def test_severity_exponent_above_one_raises():
    mg = _make_grid()
    with pytest.raises(ValidationError, match="severity_exponent must be <= 1.0"):
        WildfireGenerator(mg, severity_exponent=1.5)


def test_negative_fuel_raises():
    """fuel_availability containing negative values should raise ValueError."""
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


def test_log_not_callable_raises():
    """A non-callable log should raise TypeError."""
    mg = _make_grid()
    with pytest.raises(TypeError, match="log must be callable or None"):
        WildfireGenerator(mg, log="fires.csv")


# =============================================================================
# __init__ -- random number generator
# =============================================================================


def test_rng_none_does_not_crash():
    """Passing rng=None should initialise and run without error."""
    mg = _make_grid(fuel_val=0.5)
    wg = WildfireGenerator(mg, rng=None)
    wg.run_one_step(1)


def test_same_rng_seed_gives_same_result():
    """Two components with the same rng seed should produce identical fires."""
    _, log_a = _run_with_log(_make_grid(9, 9, fuel_val=0.8), rng=11)
    _, log_b = _run_with_log(_make_grid(9, 9, fuel_val=0.8), rng=11)
    assert log_a.records == log_b.records


def test_results_independent_of_global_numpy_seed():
    """The component must not draw from numpy's global random state."""
    np.random.seed(1)
    _, log_a = _run_with_log(_make_grid(9, 9, fuel_val=0.8), rng=11)
    np.random.seed(2)
    _, log_b = _run_with_log(_make_grid(9, 9, fuel_val=0.8), rng=11)
    assert log_a.records == log_b.records


# =============================================================================
# __init__ -- output field initialisation
# =============================================================================


def test_last_fire_time_initialised_to_neg_inf():
    """last_fire_time should be -inf at all nodes after initialisation."""
    mg = _make_grid()
    wg = WildfireGenerator(mg)
    assert np.all(
        wg._last_fire_time == -np.inf
    ), "last_fire_time should be -inf everywhere before any fire"


# =============================================================================
# __init__ -- field aliasing
# =============================================================================


def test_vegetation_field_linked_to_grid():
    """_vegetation must share memory with grid's fuel_availability field."""
    mg = _make_grid(nrows=9, ncols=9, topo_val=100, fuel_val=0.6)
    wg = WildfireGenerator(mg)
    assert np.shares_memory(wg._vegetation, mg.at_node["fuel_availability"])


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
    wg.aridity = 0.75
    assert wg.aridity == 0.75


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
    assert wg._rivers[20]
    assert not wg._rivers[21]


# =============================================================================
# FireEvent and FireEventLog
# =============================================================================


def test_fire_event_is_frozen():
    event = FireEvent(1.0, (0.0, 0.0), (6,), 1.0, 0.5, 0.5, 0.5)
    with pytest.raises(AttributeError):
        event.year = 2.0


def test_fire_event_log_starts_empty():
    log = FireEventLog()
    assert len(log) == 0
    assert log.records == []


def test_fire_event_log_collects_events_in_order():
    log = FireEventLog()
    first = FireEvent(1.0, (0.0, 0.0), (6,), 1.0, 0.5, 0.5, 0.5)
    second = FireEvent(2.0, (1.0, 1.0), (7,), 1.0, 0.5, 0.5, 0.5)
    log(first)
    log(second)
    assert len(log) == 2
    assert list(log) == [first, second]
    assert log.records == [first, second]


# =============================================================================
# Logging behaviour
# =============================================================================


def test_runs_without_log():
    """log=None (the default) should run and still burn nodes."""
    mg = _make_grid(nrows=9, ncols=9, fuel_val=0.8)
    wg = WildfireGenerator(mg, potential_fires=50, aridity=0.7)
    wg.run_one_step(1)
    assert np.any(np.isfinite(mg.at_node["last_fire_time"]))


def test_custom_callable_logger():
    """Any callable, e.g. list.append, can be used as the logger."""
    events = []
    mg = _make_grid(nrows=9, ncols=9, fuel_val=0.8)
    WildfireGenerator(
        mg, potential_fires=50, aridity=0.7, log=events.append
    ).run_one_step(1)
    assert events
    assert all(isinstance(e, FireEvent) for e in events)


def test_logging_does_not_change_physics():
    """Results must be identical with and without a logger."""
    mg_a = _make_grid(nrows=9, ncols=9, fuel_val=0.8)
    mg_b = _make_grid(nrows=9, ncols=9, fuel_val=0.8)
    WildfireGenerator(mg_a, rng=3).run_one_step(1)
    WildfireGenerator(mg_b, rng=3, log=FireEventLog()).run_one_step(1)
    testing.assert_array_equal(
        mg_a.at_node["fuel_availability"], mg_b.at_node["fuel_availability"]
    )
    testing.assert_array_equal(
        mg_a.at_node["last_fire_time"], mg_b.at_node["last_fire_time"]
    )


def test_fire_event_field_types():
    """Event fields should be plain Python types, not numpy arrays."""
    mg = _make_grid(nrows=9, ncols=9, spacing=100, fuel_val=0.8)
    _, log = _run_with_log(mg, potential_fires=50, aridity=0.7, rng=42)
    assert log.records
    for e in log:
        assert isinstance(e.xy_of_ignition, tuple)
        assert len(e.xy_of_ignition) == 2
        assert all(isinstance(v, float) for v in e.xy_of_ignition)
        assert isinstance(e.burned_nodes, tuple)
        assert all(isinstance(n, int) for n in e.burned_nodes)
        assert isinstance(e.fire_size, float)


def test_fire_sizes_positive():
    """Every fire on an all-core interior grid should have a positive size."""
    mg = _make_grid(nrows=9, ncols=9, spacing=100, fuel_val=0.8)
    _, log = _run_with_log(mg, potential_fires=50, aridity=0.7, rng=42)
    assert all(e.fire_size > 0 for e in log)


def test_severity_in_unit_range():
    mg = _make_grid(nrows=9, ncols=9, spacing=100, fuel_val=0.8)
    _, log = _run_with_log(mg, potential_fires=50, aridity=0.7, rng=42)
    assert all(0.0 <= e.severity_factor <= 1.0 for e in log)


def test_ignition_is_a_core_node():
    """Ignitions are drawn from core nodes only."""
    mg = _make_grid(nrows=9, ncols=9, spacing=100, fuel_val=0.8)
    _, log = _run_with_log(mg, potential_fires=50, aridity=0.7, rng=42)
    core_xy = {tuple(xy) for xy in mg.xy_of_node[mg.core_nodes]}
    assert all(e.xy_of_ignition in core_xy for e in log)


def test_events_from_last_step():
    """Events from the latest step all carry the current time."""
    mg = _make_grid(nrows=9, ncols=9, spacing=100, fuel_val=0.8)
    wg, log = _run_with_log(mg, n_steps=2, potential_fires=30, aridity=0.7, rng=7)
    last_step = [e for e in log if e.year == wg._current_time]
    earlier = [e for e in log if e.year != wg._current_time]
    assert all(e.year == 1.0 for e in earlier)
    assert len(last_step) + len(earlier) == len(log)


# =============================================================================
# get_regrowth_params
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
    """get_regrowth_params should return the correct (T, e) pair for each bin."""
    t, e = get_regrowth_params(aridity)
    assert t == expected_t
    assert e == expected_e


@pytest.mark.parametrize("aridity", [-0.1, 1.5])
def test_get_regrowth_params_out_of_range_raises(aridity):
    with pytest.raises(ValueError, match="outside the table range"):
        get_regrowth_params(aridity)


# =============================================================================
# get_aridity_weights
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
    """get_aridity_weights should return the correct weight pair for each bin."""
    fw, aw = get_aridity_weights(aridity)
    assert fw == exp_fw
    assert aw == exp_aw


@pytest.mark.parametrize("aridity", [0.05, 0.2, 0.4, 0.6, 0.9, 1.0])
def test_aridity_weights_sum_to_one(aridity):
    """Fuel weight + aridity weight should always sum to 1.0."""
    fw, aw = get_aridity_weights(aridity)
    assert fw + aw == pytest.approx(1.0)


@pytest.mark.parametrize("aridity", [-0.1, 1.5])
def test_get_aridity_weights_out_of_range_raises(aridity):
    with pytest.raises(ValueError, match="outside the table range"):
        get_aridity_weights(aridity)


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
    wg_low = WildfireGenerator(_make_grid(9, 9, topo_val=100), aridity=0.1, rng=0)
    wg_high = WildfireGenerator(_make_grid(9, 9, topo_val=100), aridity=0.9, rng=0)
    sev_low = [wg_low._calc_severity() for _ in range(500)]
    sev_high = [wg_high._calc_severity() for _ in range(500)]
    assert np.mean(sev_high) > np.mean(sev_low)


def test_calc_severity_at_zero_aridity_near_zero():
    """With aridity near 0, mean severity should be low."""
    wg = WildfireGenerator(_make_grid(9, 9, topo_val=100), aridity=0.01, rng=1)
    mean_sev = np.mean([wg._calc_severity() for _ in range(300)])
    assert mean_sev < 0.2, f"Mean severity at near-zero aridity was {mean_sev}"


def test_calc_severity_at_unit_aridity_near_one():
    """With aridity = 1, mean severity should be close to 1."""
    wg = WildfireGenerator(_make_grid(9, 9, topo_val=100), aridity=1.0, rng=2)
    mean_sev = np.mean([wg._calc_severity() for _ in range(300)])
    assert mean_sev > 0.9, f"Mean severity at aridity=1 was {mean_sev}"


# =============================================================================
# _regrow_vegetation
# =============================================================================


def test_regrow_vegetation_changing_dt():
    """Growth increment should match the analytic formula for each dt."""
    mg = _make_grid(nrows=5, ncols=5)
    wg = WildfireGenerator(mg)
    regrowth_time, regrowth_exponent = get_regrowth_params(wg.aridity)
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
    assert np.all(wg._vegetation >= veg_before)


def test_regrow_vegetation_capped_at_max():
    """_regrow_vegetation must not push vegetation above max_vegetation."""
    mg = _make_grid(nrows=9, ncols=9, topo_val=100, fuel_val=1.0)
    wg = WildfireGenerator(mg, aridity=0.5)
    wg._regrow_vegetation(dt=1)
    assert np.all(wg._vegetation <= wg._max_vegetation)


def test_regrow_vegetation_recovers_over_time():
    """Starting from zero, repeated regrowth should bring vegetation close to max."""
    mg = _make_grid(nrows=9, ncols=9, topo_val=100, fuel_val=0.0)
    wg = WildfireGenerator(mg, aridity=0.05)  # fastest regrowth bin
    for _ in range(500):
        wg._regrow_vegetation(dt=1)
    assert np.all(wg._vegetation > 0.9)


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
    assert np.all(wg._vegetation <= max_veg)


# =============================================================================
# _get_neighbors
# =============================================================================


def test_get_neighbors_returns_valid_nodes():
    """_get_neighbors should return only valid node indices for interior nodes."""
    mg = _make_grid(nrows=9, ncols=9, topo_val=100)
    wg = WildfireGenerator(mg)
    neighbors = wg._get_neighbors(40)  # centre of 9x9
    assert len(neighbors) > 0
    assert all(n >= 0 for n in neighbors)


def test_get_neighbors_interior_node_count():
    """A fully interior node on a 9x9 grid should have exactly 4 active D4 neighbors."""
    mg = _make_grid(nrows=9, ncols=9, topo_val=100)
    wg = WildfireGenerator(mg)
    active = [n for n in wg._get_neighbors(40) if n >= 0]
    assert len(active) == 4


# =============================================================================
# _fire
# =============================================================================


def test_fire_no_ignitions_produces_no_record():
    """When the Poisson draw returns 0, _fire should return with no events."""
    mg = _make_grid(nrows=9, ncols=9, fuel_val=0.8)
    log = FireEventLog()
    wg = WildfireGenerator(mg, potential_fires=1, rng=0, log=log)
    wg._current_time = 1
    wg._fire(dt=1e-10)  # poisson(1e-10) is almost certainly 0
    assert len(log) == 0


def test_fire_ignition_skipped_on_river_center():
    """An ignition attempt whose center node is a river should be silently skipped."""
    mg = _make_grid(nrows=9, ncols=9, topo_val=100, fuel_val=1.0)
    mg.at_node["drainage_area"][mg.core_nodes] = 1e9  # all core nodes are rivers
    log = FireEventLog()
    wg = WildfireGenerator(mg, potential_fires=200, aridity=0.9, rng=0, log=log)
    wg._current_time = 1
    wg._fire(dt=1)
    assert len(log) == 0


def test_fire_low_ignition_prob_produces_fewer_fires():
    """With very low fuel and aridity, most ignition attempts should be rejected."""
    _, log_low = _run_with_log(
        _make_grid(9, 9, spacing=100, fuel_val=0.0),
        potential_fires=50,
        aridity=0.01,
        rng=1,
    )
    _, log_high = _run_with_log(
        _make_grid(9, 9, spacing=100, fuel_val=1.0),
        potential_fires=50,
        aridity=0.9,
        rng=1,
    )
    assert len(log_high) > len(log_low)


def test_fire_reduces_vegetation():
    """Burned nodes should have lower vegetation than before the fire."""
    mg = _make_grid(nrows=9, ncols=9, topo_val=100, fuel_val=1.0)
    log = FireEventLog()
    wg = WildfireGenerator(mg, potential_fires=200, aridity=0.9, rng=7, log=log)
    wg._current_time = 1
    wg._fire(dt=1)

    assert log.records, "Expected at least one fire with high fuel/aridity/attempts"
    burned = list(log.records[0].burned_nodes)
    assert np.all(wg._vegetation[burned] < 1.0)


def test_fire_does_not_cross_rivers():
    """Fire should never burn nodes flagged as rivers."""
    mg = _make_grid(nrows=9, ncols=9, topo_val=100, fuel_val=1.0)
    river_nodes = [4, 13, 22, 31, 40, 49, 58, 67, 76]
    mg.at_node["drainage_area"][river_nodes] = 3e5  # above default threshold
    log = FireEventLog()
    wg = WildfireGenerator(mg, potential_fires=300, aridity=0.9, rng=0, log=log)
    for t in range(20):
        wg._current_time = t
        wg._fire(dt=1)

    assert log.records, "Expected fires to ignite over 20 timesteps"
    all_burned = {n for e in log for n in e.burned_nodes}
    assert all_burned.isdisjoint(river_nodes)


def test_fire_vegetation_removed_fraction_in_valid_range():
    """vegetation_removed_fraction should be in [0, 1]."""
    mg = _make_grid(nrows=9, ncols=9, spacing=100, fuel_val=0.8)
    _, log = _run_with_log(mg, potential_fires=50, aridity=0.7, rng=3)
    assert log.records
    for e in log:
        assert 0.0 <= e.vegetation_removed_fraction <= 1.0


def test_fire_vegetation_removed_fraction_equals_severity():
    """Every burned node loses the same fraction, so the mean equals severity."""
    mg = _make_grid(nrows=9, ncols=9, spacing=100, fuel_val=0.8)
    _, log = _run_with_log(mg, potential_fires=50, aridity=0.7, rng=3)
    for e in log:
        assert e.vegetation_removed_fraction == pytest.approx(e.severity_factor)


def test_fire_size_equals_cell_area_sum():
    """fire_size should equal the sum of cell areas of burned nodes, in grid units.

    Boundary nodes have cell_area = 0 so only core cells contribute.
    """
    mg = _make_grid(nrows=9, ncols=9, spacing=1000, fuel_val=1.0)
    _, log = _run_with_log(mg, potential_fires=3, aridity=0.9, rng=0)
    assert log.records
    for e in log:
        expected = mg.cell_area_at_node[list(e.burned_nodes)].sum()
        testing.assert_allclose(e.fire_size, expected, rtol=1e-10)


def test_last_fire_time_updated_at_burned_nodes():
    """last_fire_time should equal current_time at every burned node."""
    mg = _make_grid(nrows=15, ncols=15, spacing=100, fuel_val=0.8)
    wg, log = _run_with_log(mg, potential_fires=5, aridity=0.5, rng=99)
    for e in log:
        testing.assert_array_equal(
            wg._last_fire_time[list(e.burned_nodes)], wg._current_time
        )


def test_last_fire_time_unchanged_at_unburned_nodes():
    """last_fire_time should remain -inf at nodes that were never burned."""
    mg = _make_grid(nrows=15, ncols=15, spacing=100, fuel_val=0.8)
    wg, log = _run_with_log(mg, potential_fires=5, aridity=0.5, rng=99)
    burned_all = {n for e in log for n in e.burned_nodes}
    unburned = [n for n in mg.core_nodes if n not in burned_all]

    assert unburned, "Expected some core nodes to remain unburned"
    assert np.all(wg._last_fire_time[unburned] == -np.inf)


def test_fire_event_aridity_matches_component_aridity():
    """The aridity recorded in each event should match the component aridity."""
    mg = _make_grid(nrows=9, ncols=9, spacing=100, fuel_val=0.8)
    _, log = _run_with_log(mg, potential_fires=50, aridity=0.65, rng=42)
    assert log.records
    assert all(e.aridity == 0.65 for e in log)


def test_fire_revisits_already_burned_node_in_front(monkeypatch):
    """Cover the 'node in burned' branch inside the fire-spread loop."""
    mg = _make_grid(nrows=5, ncols=5, topo_val=0.0, fuel_val=1.0)
    log = FireEventLog()
    wg = WildfireGenerator(mg, potential_fires=1, aridity=0.9, log=log)

    fake_adjacency = {0: [1, 2], 1: [2], 2: [1]}
    monkeypatch.setattr(wg, "_get_neighbors", lambda node: fake_adjacency.get(node, []))
    monkeypatch.setattr(wg, "_rng", _FakeRng())

    wg._current_time = 1
    wg._fire(dt=1.0)

    assert len(log) == 1
    assert set(log.records[0].burned_nodes) == {0, 1, 2}


# =============================================================================
# run_one_step -- time tracking
# =============================================================================


def test_current_time_zero_before_first_call():
    """_current_time should be 0.0 before any call to run_one_step."""
    wg = WildfireGenerator(_make_grid())
    assert wg._current_time == 0.0


def test_run_one_step_current_time_set_on_first_call():
    """_current_time should equal dt after the first run_one_step call."""
    wg = WildfireGenerator(_make_grid())
    wg.run_one_step(5)
    assert wg._current_time == 5


def test_run_one_step_current_time_accumulates():
    """_current_time should accumulate across successive calls."""
    wg = WildfireGenerator(_make_grid())
    wg.run_one_step(2)
    assert wg._current_time == 2
    wg.run_one_step(3)
    assert wg._current_time == 5
    wg.run_one_step(7)
    assert wg._current_time == 12


# =============================================================================
# run_one_step -- fire log and vegetation
# =============================================================================


def test_fire_frequency():
    """The docstring example should produce exactly 13 fires."""
    _, log = _run_with_log(_doctest_grid())
    assert len(log) == 13


def test_fire_locations():
    """Ignition coordinates should match the seed-locked expectation."""
    _, log = _run_with_log(_doctest_grid())
    expected_locations = [
        (3.0, 2.0),
        (1.0, 1.0),
        (1.0, 2.0),
        (2.0, 2.0),
        (2.0, 1.0),
        (2.0, 2.0),
        (1.0, 1.0),
        (2.0, 3.0),
        (3.0, 1.0),
        (2.0, 3.0),
        (3.0, 3.0),
        (1.0, 1.0),
        (1.0, 3.0),
    ]
    assert [e.xy_of_ignition for e in log] == expected_locations


def test_run_one_step_updates_vegetation():
    """run_one_step should keep vegetation within [0, max_vegetation]."""
    mg = _make_grid(nrows=9, ncols=9, topo_val=100, fuel_val=0.5)
    wg = WildfireGenerator(mg, potential_fires=100, aridity=0.6, rng=42)
    wg.run_one_step(dt=1)
    assert np.all(wg._vegetation >= 0)
    assert np.all(wg._vegetation <= wg._max_vegetation)


def test_run_one_step_accumulates_events():
    """Repeated run_one_step calls should accumulate events in the log."""
    mg = _make_grid(nrows=9, ncols=9, topo_val=100, fuel_val=1.0)
    wg, log = _run_with_log(mg, n_steps=1, potential_fires=200, aridity=0.9, rng=5)
    after_one = len(log)
    for _ in range(9):
        wg.run_one_step(dt=1)
    assert len(log) > after_one


def test_run_one_step_event_year_in_simulated_range():
    """Every event's year should be within the simulated range."""
    mg = _make_grid(nrows=9, ncols=9, spacing=100, fuel_val=0.8)
    _, log = _run_with_log(mg, n_steps=5, potential_fires=30, aridity=0.7, rng=5)
    assert log.records
    for e in log:
        assert isinstance(e.year, float)
        assert 1 <= e.year <= 5

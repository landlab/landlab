import math

import numpy as np
import pytest
from requireit import ValidationError

from landlab.core.component_utils import iter_adaptive_time_steps
from landlab.core.component_utils import iter_time_steps


@pytest.mark.parametrize(
    ("duration", "dt", "expected"),
    [
        (0.0, None, []),
        (5.0, None, [5.0]),
        (10.0, 2.5, [2.5, 2.5, 2.5, 2.5]),
        (10.0, 3.0, [2.5, 2.5, 2.5, 2.5]),
        (5.0, 10.0, [5.0]),
    ],
)
def test_iter_time_steps(duration, dt, expected):
    assert list(iter_time_steps(duration, dt=dt)) == expected


@pytest.mark.parametrize("duration", [-1.0, math.inf, math.nan])
def test_iter_time_steps_rejects_invalid_duration(duration):
    with pytest.raises(ValidationError):
        list(iter_time_steps(duration))


@pytest.mark.parametrize("dt", [-1.0, 0.0, math.inf, math.nan])
def test_iter_time_steps_rejects_invalid_dt(dt):
    with pytest.raises(ValidationError):
        list(iter_time_steps(1.0, dt=dt))


def test_iter_adaptive_time_steps_caps_final_step():
    expected = [3.0, 3.0, 3.0, 1.0]
    actual = list(iter_adaptive_time_steps(10.0, calc_dt=lambda: 3.0))
    assert actual == expected


def test_iter_adaptive_time_steps_stops_on_none():
    steps = iter([2.0, 2.0, None])

    expected = [2.0, 2.0]
    actual = list(iter_adaptive_time_steps(10.0, calc_dt=lambda: next(steps)))
    assert actual == expected


def test_iter_adaptive_time_steps_accepts_infinite_step():
    assert list(iter_adaptive_time_steps(10.0, calc_dt=lambda: math.inf)) == [10.0]


def test_iter_adaptive_time_steps_does_not_calculate_step_for_zero_duration():
    def calc_dt():
        raise AssertionError("calc_dt should not be called")

    with pytest.raises(AssertionError, match="calc_dt should not be called"):
        list(iter_adaptive_time_steps(10.0, calc_dt=calc_dt))
    assert list(iter_adaptive_time_steps(0.0, calc_dt=calc_dt)) == []


def test_iter_adaptive_time_steps_stops_within_tolerance():
    step = 1.0 - 5.0e-13

    assert list(iter_adaptive_time_steps(1.0, calc_dt=lambda: step)) == [step]


def test_iter_adaptive_time_steps_honors_zero_tolerance():
    steps = list(iter_adaptive_time_steps(1.0, calc_dt=lambda: 1.0 - 5.0e-13, rtol=0.0))

    assert len(steps) == 2
    assert sum(steps) == 1.0


@pytest.mark.parametrize("max_steps", (1, np.int64(1), np.int32(1)))
def test_iter_adaptive_time_steps_accepts_int(max_steps):
    actual = list(
        iter_adaptive_time_steps(1.0, calc_dt=lambda: 1.0, max_steps=max_steps)
    )
    assert actual == [1.0]


def test_iter_adaptive_time_steps_raises_when_max_steps_exceeded():
    with pytest.raises(RuntimeError, match="unable to reach 2.0 in 1 substeps"):
        list(iter_adaptive_time_steps(2.0, calc_dt=lambda: 1.0, max_steps=1))


@pytest.mark.parametrize("max_steps", [0, -1, 1.5])
def test_iter_adaptive_time_steps_rejects_invalid_max_steps(max_steps):
    with pytest.raises(ValidationError):
        list(iter_adaptive_time_steps(1.0, calc_dt=lambda: 1.0, max_steps=max_steps))


@pytest.mark.parametrize("rtol", [-1.0, 1.0, math.inf, math.nan])
def test_iter_adaptive_time_steps_rejects_invalid_rtol(rtol):
    with pytest.raises(ValidationError):
        list(iter_adaptive_time_steps(1.0, calc_dt=lambda: 1.0, rtol=rtol))


@pytest.mark.parametrize("step", [-1.0, 0.0, math.nan])
def test_iter_adaptive_time_steps_rejects_invalid_step(step):
    with pytest.raises(RuntimeError, match="step must be positive"):
        list(iter_adaptive_time_steps(1.0, calc_dt=lambda: step))


def test_iter_adaptive_time_steps_detects_step_too_small_to_advance():
    steps = iter([1.0, 1.0e-20])

    with pytest.raises(RuntimeError, match="too small relative to the elapsed time"):
        list(iter_adaptive_time_steps(2.0, calc_dt=lambda: next(steps), rtol=0.0))

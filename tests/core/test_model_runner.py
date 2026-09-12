from dataclasses import FrozenInstanceError

import numpy as np
import pytest
from requireit import ValidationError

from landlab.core.model_runner import Clock


def test_clock_defaults():
    clock = Clock()
    assert clock.start == 0.0
    assert np.isinf(clock.stop)
    assert clock.step == 1.0
    assert np.isinf(clock.duration)


def test_clock():
    clock = Clock(start=1.0, stop=10, step=3)
    assert clock.start == 1.0
    assert clock.stop == 10
    assert clock.step == 3
    assert clock.duration == 9.0


@pytest.mark.parametrize("attr", ("start", "stop", "step"))
def test_clock_is_unchanging(attr):
    clock = Clock(start=1.0, stop=10, step=3)
    with pytest.raises(FrozenInstanceError):
        setattr(clock, attr, 999)


@pytest.mark.parametrize("step", (0.0, -1.0, np.inf, np.nan))
def test_clock_bad_step(step):
    with pytest.raises(ValueError, match="^step must"):
        Clock(step=step)


@pytest.mark.parametrize("start, stop", ((1.1, 1.0), (2.0, 2.0)))
def test_clock_start_greater_than_stop(start, stop):
    with pytest.raises(ValidationError, match="^start must be"):
        Clock(start=start, stop=stop)

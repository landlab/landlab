import numpy as np
import pytest
from numpy.testing import assert_equal

from landlab.components.overland_flow._calc import adjust_unstable_discharge


def threshold_discharge(h, dx, dt):
    return 0.25 * h * dx / dt * 0.8


@pytest.mark.parametrize("dx", (1.0, 2.0))
@pytest.mark.parametrize("dt", (1.0, 2.0))
@pytest.mark.parametrize("h", (1.0, 2.0))
def test_stable_discharge(dx, dt, h):
    eps = 1e-6

    q_max = threshold_discharge(h, dx, dt)

    q = [-q_max - eps, -q_max + eps, q_max - eps, q_max + eps]
    expected = [-q_max, -q_max + eps, q_max - eps, q_max]
    actual = np.empty_like(q)

    rtn = adjust_unstable_discharge(
        np.asarray(q),
        np.full_like(q, h),
        dx=dx,
        dt=dt,
        where=np.arange(len(q)),
        out=actual,
    )
    assert_equal(actual, expected)
    assert rtn is actual


@pytest.mark.parametrize(
    "where",
    (
        [0, 0, 0, 0],
        [1, 1, 1, 1],
        [0, 0, 0, 1],
        [1, 0, 0, 0],
    ),
)
def test_stable_discharge_where(where):
    where = np.asarray(where).astype(bool)
    dx, dt, h = 1.0, 1.0, 1.0
    eps = 1e-6

    q_max = threshold_discharge(h, dx, dt)

    q = [-q_max - eps, -q_max + eps, q_max - eps, q_max + eps]
    expected = np.asarray([-q_max, -q_max + eps, q_max - eps, q_max])
    actual = np.full_like(q, -999)

    rtn = adjust_unstable_discharge(
        np.asarray(q),
        np.full_like(q, h),
        dx=dx,
        dt=dt,
        where=np.flatnonzero(where),
        out=actual,
    )
    assert_equal(actual[where], expected[where])
    assert_equal(actual[~where], -999)
    assert rtn is actual


@pytest.mark.parametrize("dx", (1.0, 2.0))
@pytest.mark.parametrize("dt", (1.0, 2.0))
@pytest.mark.parametrize("h", (1.0, 2.0))
@pytest.mark.parametrize("n_items", (1, 5))
def test_stable_discharge_in_place(dx, dt, h, n_items):
    q_max = threshold_discharge(h, dx, dt)

    h = np.full(n_items, h, dtype=float)

    q = np.full_like(h, q_max * (1 + 1e-6))
    out = q
    q_initial = q.copy()
    expected = np.full_like(h, q_max)

    rtn = adjust_unstable_discharge(
        q, h, dx=dx, dt=dt, where=np.arange(n_items), out=out
    )

    assert rtn is out
    assert_equal(out, expected)
    assert_equal(q, expected)
    assert np.all(q != q_initial)

import numpy as np
import pytest
from numpy.testing import assert_equal

from landlab.components.overland_flow._calc import adjust_supercritical_discharge


def supercritical_reference(h, g, froude):
    return froude * h * np.sqrt(g * h)


@pytest.mark.parametrize("h", (0.0, 1.0, 2.0))
@pytest.mark.parametrize("g", (0.0, 9.81))
@pytest.mark.parametrize("froude", (0.0, 1.0, 2.0))
def test_supercritical(h, g, froude):
    q_max = supercritical_reference(h, g, froude)

    q = np.asarray([-2 * q_max, -q_max, 0.0, q_max, 2 * q_max], dtype=float)
    h = np.full_like(q, h)
    where = np.arange(len(q), dtype=int)
    actual = np.empty_like(q)

    expected = np.clip(q, -q_max, q_max)

    rtn = adjust_supercritical_discharge(
        q, h, g=g, froude=froude, where=where, out=actual
    )
    assert rtn is actual
    assert_equal(actual, expected)
    assert_equal(np.sign(actual), np.sign(q))


@pytest.mark.parametrize("h", (1.0,))
@pytest.mark.parametrize("g", (9.81,))
@pytest.mark.parametrize("froude", (1.0,))
def test_supercritical_in_place(h, g, froude):
    q_max = supercritical_reference(h, g, froude)

    q = np.asarray([-2 * q_max, -q_max, 0.0, q_max, 2 * q_max], dtype=float)
    h = np.full_like(q, h)
    where = np.arange(len(q), dtype=int)
    actual = q
    expected = np.clip(q, -q_max, q_max)
    rtn = adjust_supercritical_discharge(
        q, h, g=g, froude=froude, where=where, out=actual
    )
    assert rtn is actual
    assert_equal(actual, expected)
    assert_equal(actual, q)


@pytest.mark.parametrize(
    "where",
    (
        [0, 0, 0, 0, 0],
        [1, 1, 1, 1, 1],
        [0, 1, 0, 1, 0],
        [0, 0, 0, 0, 1],
    ),
)
@pytest.mark.parametrize("h", (1.0,))
@pytest.mark.parametrize("g", (9.81,))
@pytest.mark.parametrize("froude", (1.0,))
def test_supercritical_where(where, h, g, froude):
    q_max = supercritical_reference(h, g, froude)

    q = np.asarray([-2 * q_max, -q_max, 0.0, q_max, 2 * q_max], dtype=float)
    h = np.full_like(q, h)
    where = np.asarray(where).astype(bool)
    actual = np.full_like(q, -999)
    expected = np.clip(q, -q_max, q_max)

    rtn = adjust_supercritical_discharge(
        q, h, g=g, froude=froude, where=np.flatnonzero(where), out=actual
    )
    assert rtn is actual
    assert np.all(actual[~where] == -999)
    assert_equal(actual[where], expected[where])

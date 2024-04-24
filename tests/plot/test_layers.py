import hypothesis.extra.numpy as hynp
import pytest
from hypothesis import given
from hypothesis.strategies import floats
from hypothesis.strategies import integers
from numpy.testing import assert_array_equal

from landlab.plot.layers import _insert_shorelines
from landlab.plot.layers import _interp_zero_crossings
from landlab.plot.layers import _search_zero_crossings


@given(
    y=hynp.arrays(
        dtype=hynp.floating_dtypes(),
        shape=8,
        elements=floats(0, 1, width=16, exclude_min=False),
    )
)
def test_search_zero_crossings_no_shore_floats(y):
    indices = _search_zero_crossings(y)
    assert len(indices) == 0

    indices = _search_zero_crossings(y * -1)
    assert len(indices) == 0


@given(
    y=hynp.arrays(
        dtype=int,
        shape=16,
        elements=integers(0, 128),
    )
)
def test_search_zero_crossings_no_shore_ints(y):
    indices = _search_zero_crossings(y)
    assert len(indices) == 0

    indices = _search_zero_crossings(y * -1)
    assert len(indices) == 0


def test_search_zero_crossings_all_zeros():
    assert len(_search_zero_crossings([0] * 128)) == 0


@pytest.mark.parametrize(
    "y,shoreline",
    [
        ([1, -1, -1], [0]),
        ([1, 0, 1], []),
        ([1, 0, -1], []),
        ([1, 0, 0, -1], []),
        ([-2, 0, 0, 1], []),
        ([1, -1, 1, -1], [0, 1, 2]),
    ],
)
def test_search_zero_crossings(y, shoreline):
    assert_array_equal(_search_zero_crossings(y), shoreline)


def test_interp_zero_crossings():
    assert_array_equal(_interp_zero_crossings([0, 1, 2], [1, -1, -1], [0]), [0.5])
    assert_array_equal(
        _interp_zero_crossings([0, 1, 2, 3], [1, -1, -1, 4], [0, 2]), [0.5, 2.2]
    )


def test_interp_zero_crossings_on_end_point():
    assert_array_equal(_interp_zero_crossings([1.0, 3.0, 5.0], [1, 0, -1], [1]), [3.0])
    assert_array_equal(_interp_zero_crossings([1.0, 3.0, 5.0], [1, 0, -1], [0]), [3.0])


@pytest.mark.parametrize("shoreline", [0, 1])
def test_interp_zero_crossings_bounds_error(shoreline):
    with pytest.raises(ValueError):
        _interp_zero_crossings([1.0, 3.0, 5.0], [1, 2, 3], [shoreline])


def test_interp_zero_crossings_last_element():
    assert_array_equal(_interp_zero_crossings([1.0, 3.0, 5.0], [1, 1, 0], [2]), [5.0])


def test_interp_shoreline_all_shoreline():
    assert_array_equal(_interp_zero_crossings([1.0, 3.0, 5.0], [0, 0, 0], [0]), [3.0])
    assert_array_equal(_interp_zero_crossings([1.0, 3.0, 5.0], [0, 0, 0], [1]), [5.0])
    assert_array_equal(
        _interp_zero_crossings([1.0, 3.0, 5.0], [0, 0, 0], [0, 1]), [3.0, 5.0]
    )


def test_insert_shorelines():
    x, y = _insert_shorelines([0, 1, 2], [1, -1, -1])
    assert_array_equal(x, [0, 0.5, 1, 2])
    assert_array_equal(y, [1, 0.0, -1, -1])


def test_insert_shorelines_already_there():
    x0, y0 = [0, 1, 2], [1, 0, -1]
    x, y = _insert_shorelines([0, 1, 2], [1, 0, -1])
    assert_array_equal(x, x0)
    assert_array_equal(y, y0)


def test_insert_shorelines_multiple():
    x, y = _insert_shorelines([0, 1, 2, 3, 4, 5], [1, -1, -1, 3, 4, 5])
    assert_array_equal(x, [0, 0.5, 1, 2, 2.25, 3, 4, 5])
    assert_array_equal(y, [1, 0.0, -1, -1, 0.0, 3, 4, 5])


@given(
    y0=hynp.arrays(
        dtype=hynp.floating_dtypes(),
        shape=8,
        elements=floats(0, 1, width=16, exclude_min=False),
    )
)
def test_insert_shorelines_none(y0):
    x0 = range(len(y0))
    x, y = _insert_shorelines(x0, y0)
    assert_array_equal(x, x0)
    assert_array_equal(y, y0)

    y0 *= -1
    x, y = _insert_shorelines(x0, y0)
    assert_array_equal(x, x0)
    assert_array_equal(y, y0)

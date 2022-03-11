from landlab.plot.layers import _interp_shorelines, _search_shorelines
from hypothesis.strategies import integers, floats
from hypothesis import given
from numpy.testing import assert_array_equal
import hypothesis.extra.numpy as hynp
import pytest


@given(
    y=hynp.arrays(
        dtype=hynp.floating_dtypes(),
        shape=8,
        elements=floats(0, 1, width=16),
    )
)
def test_search_shorelines_no_shore_floats(y):
    indices = _search_shorelines(y)
    assert len(indices) == 0

    indices = _search_shorelines(y * -1)
    assert len(indices) == 0


@given(
    y=hynp.arrays(
        dtype=int,
        shape=16,
        elements=integers(0, 128),
    )
)
def test_search_shorelines_no_shore_ints(y):
    indices = _search_shorelines(y)
    assert len(indices) == 0

    indices = _search_shorelines(y * -1)
    assert len(indices) == 0


def test_search_shorelines_all_zeros():
    assert len(_search_shorelines([0] * 128)) == 0


@pytest.mark.parametrize(
    "y,shoreline",
    [
        ([1, -1, -1], [0]),
        ([1, 0, 1], []),
        ([1, 0, -1], [1]),
        ([1, 0, 0, -1], [2]),
        ([-2, 0, 0, 1], [2]),
    ],
)
def test_search_shorelines(y, shoreline):
    assert_array_equal(_search_shorelines(y), shoreline)


def test_interp_shorelines():
    assert_array_equal(_interp_shorelines([0, 1, 2], [1, -1, -1], [0]), [0.5])
    assert_array_equal(
        _interp_shorelines([0, 1, 2, 3], [1, -1, -1, 4], [0, 2]), [0.5, 2.2]
    )


def test_interp_shorelines_on_end_point():
    assert_array_equal(_interp_shorelines([1.0, 3.0, 5.0], [1, 0, -1], [1]), [3.0])
    assert_array_equal(_interp_shorelines([1.0, 3.0, 5.0], [1, 0, -1], [0]), [3.0])


@pytest.mark.parametrize("shoreline", [0, 1])
def test_interp_shorelines_bounds_error(shoreline):
    with pytest.raises(ValueError):
        _interp_shorelines([1.0, 3.0, 5.0], [1, 2, 3], [shoreline])


def test_interp_shorelines_last_element():
    with pytest.raises(ValueError):
        _interp_shorelines([1.0, 3.0, 5.0], [1, 1, 0], [2])


def test_interp_shoreline_all_shoreline():
    assert_array_equal(_interp_shorelines([1.0, 3.0, 5.0], [0, 0, 0], [1]), [5.0])

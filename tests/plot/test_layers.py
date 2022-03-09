from landlab.plot.layers import _search_shorelines
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

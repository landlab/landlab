import numpy as np
import pytest
from numpy.testing import assert_array_equal

from landlab.utils.geometry.planar import calc_cumulative_path_length
from landlab.utils.geometry.planar import calc_path_length
from landlab.utils.geometry.planar import calc_path_segment_lengths


@pytest.mark.parametrize(
    "points",
    (
        [[0.0, 0.0], [1.0, 0.0]],
        [[0.0, 0.0], [0.0, -1.0]],
        [[1.0, 2.0], [1.5, 2.0], [1.5, 2.5]],
        [[1.0, 0.0, 2.0], [1.5, 0.0, 2.0], [1.5, 0.0, 2.5]],
    ),
)
def test_path_length(points):
    assert calc_path_length(points) == 1.0
    assert_array_equal(calc_cumulative_path_length(points)[[0, -1]], (0.0, 1.0))
    assert len(calc_cumulative_path_length(points)) == len(points)


def test_path_with_one_point():
    points = [[1.0, 2.0]]
    assert_array_equal(calc_path_segment_lengths(points), [])
    assert_array_equal(calc_cumulative_path_length(points), [0.0])
    assert calc_path_length(points) == 0.0


@pytest.mark.parametrize("points", (np.empty((0, 2)), np.empty((0, 3)), []))
def test_path_with_no_points(points):
    assert_array_equal(calc_path_segment_lengths(points), [])
    assert_array_equal(calc_cumulative_path_length(points), [])
    assert calc_path_length(points) == 0.0


def test_cumulative_path_length():
    assert_array_equal(
        calc_cumulative_path_length([[0.0, 0.0], [1.0, 0.0], [1.0, 1.0], [2.0, 1.0]]),
        [0.0, 1.0, 2.0, 3.0],
    )


def test_path_segment_lengths():
    assert_array_equal(
        calc_path_segment_lengths([[0.0, 0.0], [1.0, 0.0], [1.0, 2.0], [4.0, 2.0]]),
        [1.0, 2.0, 3.0],
    )


@pytest.mark.parametrize("ndim", (1, 2, 3, 4))
def test_wrong_dimensionality(ndim):
    points = np.zeros((1,) * ndim)
    if ndim == 2:
        assert_array_equal(calc_path_segment_lengths(points), [])
    else:
        with pytest.raises(ValueError) as excinfo:
            calc_path_segment_lengths(points)
        assert str(excinfo.value).startswith("Expected a 2D array")

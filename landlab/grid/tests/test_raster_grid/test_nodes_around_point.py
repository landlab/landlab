import numpy as np
from numpy.testing import assert_array_equal

from landlab import RasterModelGrid


def test_lower_left_cell():
    rmg = RasterModelGrid((3, 3))
    assert_array_equal(rmg.nodes_around_point(0.1, 0.1), np.array([0, 3, 4, 1]))


def test_upper_left_cell():
    rmg = RasterModelGrid((3, 3))
    assert_array_equal(rmg.nodes_around_point(0.1, 0.9), np.array([0, 3, 4, 1]))


def test_upper_right_cell():
    rmg = RasterModelGrid((3, 3))
    assert_array_equal(rmg.nodes_around_point(0.9, 0.9), np.array([0, 3, 4, 1]))


def test_lower_right_cell():
    rmg = RasterModelGrid((3, 3))
    assert_array_equal(rmg.nodes_around_point(0.9, 0.1), np.array([0, 3, 4, 1]))


def test_on_left_edge():
    rmg = RasterModelGrid((3, 3))
    assert_array_equal(rmg.nodes_around_point(0.0, 0.5), np.array([0, 3, 4, 1]))


def test_on_bottom_edge():
    rmg = RasterModelGrid((3, 3))
    assert_array_equal(rmg.nodes_around_point(0.5, 0.0), np.array([0, 3, 4, 1]))


def test_on_right_edge():
    rmg = RasterModelGrid((3, 3))
    assert_array_equal(rmg.nodes_around_point(1.0, 0.5), np.array([1, 4, 5, 2]))


def test_on_top_edge():
    rmg = RasterModelGrid((3, 3))
    assert_array_equal(rmg.nodes_around_point(0.5, 1.0), np.array([3, 6, 7, 4]))


def test_on_lower_left_corner():
    rmg = RasterModelGrid((3, 3))
    assert_array_equal(rmg.nodes_around_point(0.0, 0.0), np.array([0, 3, 4, 1]))


def test_on_upper_right_corner():
    rmg = RasterModelGrid((3, 3))
    assert_array_equal(rmg.nodes_around_point(1.0, 1.0), np.array([4, 7, 8, 5]))


def test_on_upper_left_corner():
    rmg = RasterModelGrid((3, 3))
    assert_array_equal(rmg.nodes_around_point(0.0, 1.0), np.array([3, 6, 7, 4]))


def test_on_lower_right_corner():
    rmg = RasterModelGrid((3, 3))
    assert_array_equal(rmg.nodes_around_point(1.0, 0.0), np.array([1, 4, 5, 2]))

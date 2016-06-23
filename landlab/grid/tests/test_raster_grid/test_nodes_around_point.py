import numpy as np
from numpy.testing import assert_array_equal
from nose.tools import with_setup


def setup_grid():
    from landlab import RasterModelGrid
    globals().update({
        'rmg': RasterModelGrid(3, 3)
    })


@with_setup(setup_grid)
def test_lower_left_cell():
    assert_array_equal(rmg.nodes_around_point(.1, .1),
                       np.array([0, 3, 4, 1]))


@with_setup(setup_grid)
def test_upper_left_cell():
    assert_array_equal(rmg.nodes_around_point(.1, .9),
                       np.array([0, 3, 4, 1]))


@with_setup(setup_grid)
def test_upper_right_cell():
    assert_array_equal(rmg.nodes_around_point(.9, .9),
                       np.array([0, 3, 4, 1]))


@with_setup(setup_grid)
def test_lower_right_cell():
    assert_array_equal(rmg.nodes_around_point(.9, .1),
                       np.array([0, 3, 4, 1]))


@with_setup(setup_grid)
def test_on_left_edge():
    assert_array_equal(rmg.nodes_around_point(0., .5),
                       np.array([0, 3, 4, 1]))


@with_setup(setup_grid)
def test_on_bottom_edge():
    assert_array_equal(rmg.nodes_around_point(.5, 0.),
                       np.array([0, 3, 4, 1]))


@with_setup(setup_grid)
def test_on_right_edge():
    assert_array_equal(rmg.nodes_around_point(1., .5),
                       np.array([1, 4, 5, 2]))


@with_setup(setup_grid)
def test_on_top_edge():
    assert_array_equal(rmg.nodes_around_point(.5, 1.),
                       np.array([3, 6, 7, 4]))


@with_setup(setup_grid)
def test_on_lower_left_corner():
    assert_array_equal(rmg.nodes_around_point(0., 0.),
                       np.array([0, 3, 4, 1]))


@with_setup(setup_grid)
def test_on_upper_right_corner():
    assert_array_equal(rmg.nodes_around_point(1., 1.),
                       np.array([4, 7, 8, 5]))


@with_setup(setup_grid)
def test_on_upper_left_corner():
    assert_array_equal(rmg.nodes_around_point(0., 1.),
                       np.array([3, 6, 7, 4]))


@with_setup(setup_grid)
def test_on_lower_right_corner():
    assert_array_equal(rmg.nodes_around_point(1., 0.),
                       np.array([1, 4, 5, 2]))

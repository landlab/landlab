import numpy as np
from numpy.testing import assert_array_equal, assert_raises
from nose.tools import assert_equal, assert_tuple_equal

from landlab import RasterModelGrid


def test_nodes_at_left_edge():
    grid = RasterModelGrid((3, 4))
    assert_array_equal(grid.nodes_at_left_edge,
                       np.array([0, 4, 8], dtype=np.int))

    vals_at_node = np.array([ 0,  1,  2,  3,
                              4,  5,  6,  7,
                              8,  9, 10, 11])
    assert_array_equal(vals_at_node[grid.nodes_at_left_edge],
                       np.array([0, 4, 8]))


def test_nodes_at_right_edge():
    grid = RasterModelGrid((3, 4))
    assert_array_equal(grid.nodes_at_right_edge,
                       np.array([3, 7, 11], dtype=np.int))


def test_nodes_at_top_edge():
    grid = RasterModelGrid((3, 4))
    assert_array_equal(grid.nodes_at_top_edge,
                       np.array([8, 9, 10, 11], dtype=np.int))


def test_nodes_at_bottom_edge():
    grid = RasterModelGrid((3, 4))
    assert_array_equal(grid.nodes_at_bottom_edge,
                       np.array([0, 1, 2, 3], dtype=np.int))


def test_nodes_at_edge():
    grid = RasterModelGrid((3, 4))
    for edge in ('right', 'top', 'left', 'bottom'):
        assert_array_equal(grid.nodes_at_edge(edge),
                           getattr(grid, 'nodes_at_{0}_edge'.format(edge)))
    with assert_raises(ValueError):
        grid.nodes_at_edge('not-an-edge')


def test_grid_shape():
    grid = RasterModelGrid((3, 4))
    assert_tuple_equal(grid.shape, (3, 4))


def test_node_rows():
    grid = RasterModelGrid((4, 5))
    assert_equal(grid.number_of_node_rows, 4)


def test_node_columns():
    grid = RasterModelGrid((4, 5))
    assert_equal(grid.number_of_node_columns, 5)


def test_nodes_at_patch():
    grid = RasterModelGrid((3, 3))
    assert_array_equal(grid.nodes_at_patch,
                       np.array([[4, 3, 0, 1], [5, 4, 1, 2],
                                 [7, 6, 3, 4], [8, 7, 4, 5]], dtype=int))
    with assert_raises(ValueError):
        grid.nodes_at_patch[0, 0] = 42

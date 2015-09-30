import numpy as np
from numpy.testing import assert_array_equal
from nose.tools import with_setup

from landlab import RasterModelGrid
from landlab.grid.base import BAD_INDEX_VALUE


def test_all_neighbors():
    rmg = RasterModelGrid(5, 4)
    X = BAD_INDEX_VALUE
    expected = np.array([
        [X, X, X, X], [X, 5, X, X], [X, 6, X, X], [X, X, X, X],
        [5, X, X, X], [6, 9, 4, 1], [7, 10, 5, 2], [X, X, 6, X],
        [9, X, X, X], [10, 13, 8, 5], [11, 14, 9, 6], [X, X, 10, X],
        [13, X, X, X], [14, 17, 12, 9], [15, 18, 13, 10], [X, X, 14, X],
        [X, X, X, X], [X, X, X, 13], [X, X, X, 14], [X, X, X, X],
    ])
    assert_array_equal(rmg.get_neighbor_list(), expected)


def test_neighbor_list_with_scalar_arg():
    X = BAD_INDEX_VALUE
    rmg = RasterModelGrid(5, 4)

    assert_array_equal(rmg.get_neighbor_list(6), np.array([7, 10, 5, 2]))
    assert_array_equal(rmg.get_neighbor_list(-1), np.array([X, X, X, X]))
    assert_array_equal(rmg.get_neighbor_list(-2), np.array([X, X, X, 14]))


def test_neighbor_list_with_array_arg():
    X = BAD_INDEX_VALUE
    rmg = RasterModelGrid(5, 4)
    assert_array_equal(rmg.get_neighbor_list([6, -1]),
                       np.array([[7, 10, 5, 2], [X, X, X, X]]))


def test_neighbor_list_boundary():
    """
    All of the neighbor IDs for a boundary cell are -1.
    """
    X = BAD_INDEX_VALUE
    rmg = RasterModelGrid(5, 4)
    import landlab.utils.structured_grid as sgrid
    rmg.set_closed_nodes([0, 1, 2, 3, 4, 7, 8, 11, 12, 15, 16, 17, 18, 19])

    for node_id in sgrid.boundary_iter(rmg.shape):
        assert_array_equal(rmg.get_neighbor_list(node_id),
                           np.array([X, X, X, X]))

import numpy as np
from numpy.testing import assert_array_equal
from nose.tools import with_setup

from landlab import RasterModelGrid
from landlab import FIXED_VALUE_BOUNDARY as FV
from landlab import CLOSED_BOUNDARY as CB


def test_get_status():
    """Test getting the node status array."""
    grid = RasterModelGrid((4, 5))
    assert_array_equal(grid.status_at_node,
                       [FV, FV, FV, FV, FV,
                        FV,  0,  0,  0, FV,
                        FV,  0,  0,  0, FV,
                        FV, FV, FV, FV, FV])


def test_set_status_with_scalar():
    """Test setting status with a scalar."""
    grid = RasterModelGrid((4, 5))
    grid.status_at_node[6] = 2

    assert_array_equal(grid.status_at_node,
                       [FV, FV, FV, FV, FV,
                        FV,  2,  0,  0, FV,
                        FV,  0,  0,  0, FV,
                        FV, FV, FV, FV, FV])


def test_set_status_with_slice():
    """Test setting status with an array slice."""
    grid = RasterModelGrid((4, 5))
    grid.status_at_node[6: 9] = 2

    assert_array_equal(grid.status_at_node,
                       [FV, FV, FV, FV, FV,
                        FV,  2,  2,  2, FV,
                        FV,  0,  0,  0, FV,
                        FV, FV, FV, FV, FV])


def test_set_status_with_array_bool():
    """Test setting node status with boolean array."""
    grid = RasterModelGrid((4, 5))
    inds = np.full((20, ), False, dtype=bool)
    inds[6] = True
    inds[7] = True
    inds[13] = True

    grid.status_at_node[inds] = 2
    assert_array_equal(grid.status_at_node,
                       [FV, FV, FV, FV, FV,
                        FV,  2,  2,  0, FV,
                        FV,  0,  0,  2, FV,
                        FV, FV, FV, FV, FV])


def test_set_with_itemset():
    grid = RasterModelGrid((4, 5))
    grid.status_at_node.itemset(7, 2)

    assert_array_equal(grid.status_at_node,
                       [FV, FV, FV, FV, FV,
                        FV,  0,  2,  0, FV,
                        FV,  0,  0,  0, FV,
                        FV, FV, FV, FV, FV])


def test_set_status_with_array():
    """Test that active links are reset after changing the node status."""
    grid = RasterModelGrid((4, 5))

    assert_array_equal(grid.active_links,
                      [ 5,  6,  7,  9, 10, 11, 12, 14, 15,
                       16, 18, 19, 20, 21, 23, 24, 25])

    grid.status_at_node[: 5] = CB
    assert_array_equal(grid.status_at_node,
                       [CB, CB, CB, CB, CB,
                        FV,  0,  0,  0, FV,
                        FV,  0,  0,  0, FV,
                        FV, FV, FV, FV, FV])

    assert_array_equal(grid.active_links,
                      [ 9, 10, 11, 12, 14, 15, 16, 18, 19, 20, 21, 23, 24, 25])

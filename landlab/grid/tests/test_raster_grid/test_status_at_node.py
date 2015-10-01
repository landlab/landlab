import numpy as np
from numpy.testing import assert_array_equal
from nose.tools import with_setup

from landlab import RasterModelGrid
from landlab import FIXED_VALUE_BOUNDARY as FV


def test_get_status():
    """Test getting the node status array."""
    grid = RasterModelGrid((4, 5))
    assert_array_equal(grid.status_at_node,
                       [FV, FV, FV, FV, FV,
                        FV,  0,  0,  0, FV,
                        FV,  0,  0,  0, FV,
                        FV, FV, FV, FV, FV])


def test_set_status_with_scalar():
    grid = RasterModelGrid((4, 5))
    grid.status_at_node[6] = 2

    assert_array_equal(grid.status_at_node,
                       [FV, FV, FV, FV, FV,
                        FV,  2,  0,  0, FV,
                        FV,  0,  0,  0, FV,
                        FV, FV, FV, FV, FV])


def test_set_status_with_slice():
    grid = RasterModelGrid((4, 5))
    grid.status_at_node[6: 9] = 2

    assert_array_equal(grid.status_at_node,
                       [FV, FV, FV, FV, FV,
                        FV,  2,  2,  2, FV,
                        FV,  0,  0,  0, FV,
                        FV, FV, FV, FV, FV])


def test_set_status_with_array():
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

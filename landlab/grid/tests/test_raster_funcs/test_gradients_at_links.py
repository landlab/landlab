import numpy as np
from numpy.testing import assert_array_equal
from nose import with_setup
try:
    from nose.tools import assert_is
except ImportError:
    from landlab.testing.tools import assert_is

from landlab import RasterModelGrid


_GRID = None
_VALUES_AT_NODES = None


def setup_unit_grid():
    """Set up a test grid with unit spacing."""
    global _GRID, _VALUES_AT_NODES
    _GRID = RasterModelGrid(4, 5)
    _VALUES_AT_NODES = np.arange(20.)


def setup_non_unit_grid():
    """Set up a test grid with non-unit spacing."""
    global _GRID, _VALUES_AT_NODES
    _GRID = RasterModelGrid(4, 5, 2.)
    _VALUES_AT_NODES = np.arange(20.)


def setup_non_square_grid():
    """Set up a test grid with different x and y spacing."""
    global _GRID, _VALUES_AT_NODES
    _GRID = RasterModelGrid((4, 5), spacing=(5, 2))
    _VALUES_AT_NODES = np.arange(20.)


@with_setup(setup_unit_grid)
def test_unit_spacing():
    """Test on a grid with unit spacing."""
    grads = _GRID.calc_grad_at_link(_VALUES_AT_NODES)
    assert_array_equal(
        grads,
        np.array([1, 1, 1, 1, 
                  5, 5, 5, 5, 5,
                  1, 1, 1, 1,
                  5, 5, 5, 5, 5,
                  1, 1, 1, 1,
                  5, 5, 5, 5, 5,
                  1, 1, 1, 1],
                 dtype=float))
    diffs = _GRID.calculate_diff_at_links(_VALUES_AT_NODES)
    assert_array_equal(grads, diffs)


@with_setup(setup_non_square_grid)
def test_non_unit_spacing():
    """Test on a grid with non-unit spacing."""
    grads = _GRID.calc_grad_at_link(_VALUES_AT_NODES)
    assert_array_equal(
        grads,
        np.array(
           [0.5, 0.5,  0.5,  0.5,
            1.0, 1.0,  1. ,  1. ,  1. ,  
            0.5,  0.5, 0.5,  0.5,
            1. ,  1. , 1. ,  1. ,  1. ,
            0.5,  0.5,  0.5,  0.5,
            1. ,  1. ,  1. ,  1. , 1. ,
            0.5,  0.5,  0.5,  0.5],
            dtype=float))
    diffs = _GRID.calc_diff_at_link(_VALUES_AT_NODES)
    assert_array_equal(diffs,
        np.array([1, 1, 1, 1,
                  5, 5, 5, 5, 5,
                  1, 1, 1, 1,
                  5, 5, 5, 5, 5,
                  1, 1, 1, 1,
                  5, 5, 5, 5, 5,
                  1, 1, 1, 1],
                 dtype=float))


@with_setup(setup_unit_grid)
def test_out_array():
    """Test using the out keyword."""
    grads = np.empty(31)
    rtn_grads = _GRID.calc_grad_at_link(_VALUES_AT_NODES, out=grads)
    assert_array_equal(
        grads,
        np.array([1, 1, 1, 1, 5, 5, 5, 5, 5, 1, 1, 1, 1, 5, 5,
                  5, 5, 5, 1, 1, 1, 1, 5, 5, 5, 5, 5, 1, 1, 1, 1],
                 dtype=float))
    assert_is(rtn_grads, grads)


@with_setup(setup_unit_grid)
def test_diff_out_array():
    """Test that return array is the out keyword."""
    diff = np.empty(31)
    rtn_diff = _GRID.calculate_diff_at_links(_VALUES_AT_NODES, out=diff)
    assert_array_equal(
        diff,
        np.array([1, 1, 1, 1, 5, 5, 5, 5, 5, 1, 1, 1, 1, 5, 5,
                  5, 5, 5, 1, 1, 1, 1, 5, 5, 5, 5, 5, 1, 1, 1, 1],
                 dtype=float))
    assert_is(rtn_diff, diff)

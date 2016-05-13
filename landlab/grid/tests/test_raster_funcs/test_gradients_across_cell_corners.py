import numpy as np
from numpy.testing import assert_array_equal
from nose import with_setup
from nose.tools import assert_equal
try:
    from nose.tools import assert_is
except ImportError:
    from landlab.testing.tools import assert_is

from landlab.grid.raster_gradients import (
    calc_grad_across_cell_corners)


def setup_unit_grid():
    """Set up a test grid with unit spacing."""
    from landlab import RasterModelGrid
    globals().update({
        'rmg': RasterModelGrid(4, 5),
        'values_at_nodes':  np.arange(20.),
    })


def setup_grid():
    """Set up a test grid."""
    from landlab import RasterModelGrid
    globals().update({
        'rmg': RasterModelGrid(4, 5, 2.),
        'values_at_nodes':  np.arange(20.),
    })


@with_setup(setup_unit_grid)
def test_scalar_arg():
    """Test using a scalar for cell arg."""
    grads = calc_grad_across_cell_corners(
        rmg, values_at_nodes, 0)
    assert_array_equal(grads, np.array([[6., 4., -6., -4.]]) / np.sqrt(2.))


@with_setup(setup_unit_grid)
def test_iterable():
    """Test using an iterable for cell arg."""
    grads = rmg.calc_grad_across_cell_corners(values_at_nodes, [0, 4])
    assert_array_equal(grads, np.array([[6., 4., -6., -4.],
                                        [6., 4., -6., -4.]]) / np.sqrt(2.))


@with_setup(setup_unit_grid)
def test_with_no_cell_id_arg():
    """Test without using an arg for cell id."""
    values = np.array([0, 1,  3, 6, 10,
                       0, 1,  3, 6, 10,
                       0, 1,  3, 5, 10,
                       0, 1, -3, 6, 10], dtype=float)
    grads = rmg.calc_grad_across_cell_corners(values)

    assert_array_equal(grads, (1. / np.sqrt(2.)) * np.array([
        [2., -1., -1., 2.], [2., -2., -2., 3.], [4., -3., -3., 4.],
        [-4., -1., -1., 2.], [3., -2., -2., 3.], [5., -8., -2., 5.]]))


@with_setup(setup_unit_grid)
def test_with_out_keyword():
    """Test with out keyword."""
    out = np.empty((1, 4))
    rtn = rmg.calc_grad_across_cell_corners(
        values_at_nodes, 5, out=out)
    assert_is(rtn, out)
    assert_array_equal(out, np.array([[6., 4., -6., -4.]]) / np.sqrt(2))

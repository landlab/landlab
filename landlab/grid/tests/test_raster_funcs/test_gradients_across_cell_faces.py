import numpy as np
from numpy.testing import assert_array_equal
from nose import with_setup
from nose.tools import assert_equal
try:
    from nose.tools import assert_is
except ImportError:
    from landlab.testing.tools import assert_is

from landlab.grid.raster_gradients import calc_grad_across_cell_faces


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
    """Test with a scalar arg for faces."""
    grads = calc_grad_across_cell_faces(
        rmg, values_at_nodes, 0)
    assert_array_equal(grads, np.array([[1., 5., -1., -5.]]))


@with_setup(setup_unit_grid)
def test_iterable():
    """Test with an iterable arg for faces."""
    grads = rmg.calc_grad_across_cell_faces(values_at_nodes, [0, 4])
    assert_array_equal(grads, np.array([[1., 5., -1., -5.],
                                        [1., 5., -1., -5.]]))


@with_setup(setup_unit_grid)
def test_with_no_cell_id_arg():
    """Test without an arg for faces."""
    values = np.array([0, 1,  3, 6, 10,
                       0, 1,  3, 6, 10,
                       0, 1,  3, 5, 10,
                       0, 1, -3, 6, 10], dtype=float)
    grads = rmg.calc_grad_across_cell_faces(values)

    assert_array_equal(grads, np.array([
        [2., 0., -1., 0.], [3.,  0., -2., 0.], [4., -1., -3., 0.],
        [2., 0., -1., 0.], [2., -6., -2., 0.], [5.,  1., -2., 1.]]))


@with_setup(setup_unit_grid)
def test_with_out_keyword():
    """Test using the out keyword."""
    out = np.empty((1, 4))
    rtn = rmg.calc_grad_across_cell_faces(values_at_nodes, 5, out=out)
    assert_is(rtn, out)
    assert_array_equal(out, np.array([[1., 5., -1., -5.]]))

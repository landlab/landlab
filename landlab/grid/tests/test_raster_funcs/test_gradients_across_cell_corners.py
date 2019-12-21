import numpy as np
from numpy.testing import assert_array_equal

from landlab import RasterModelGrid
from landlab.grid.raster_gradients import calc_grad_across_cell_corners


def test_scalar_arg():
    """Test using a scalar for cell arg."""
    rmg = RasterModelGrid((4, 5))
    values_at_nodes = np.arange(20.0)
    grads = calc_grad_across_cell_corners(rmg, values_at_nodes, 0)
    assert_array_equal(grads, np.array([[6.0, 4.0, -6.0, -4.0]]) / np.sqrt(2.0))


def test_iterable():
    """Test using an iterable for cell arg."""
    rmg = RasterModelGrid((4, 5))
    values_at_nodes = np.arange(20.0)
    grads = rmg.calc_grad_across_cell_corners(values_at_nodes, [0, 4])
    assert_array_equal(
        grads, np.array([[6.0, 4.0, -6.0, -4.0], [6.0, 4.0, -6.0, -4.0]]) / np.sqrt(2.0)
    )


def test_with_no_cell_id_arg():
    """Test without using an arg for cell id."""
    rmg = RasterModelGrid((4, 5))
    values = np.array(
        [0, 1, 3, 6, 10, 0, 1, 3, 6, 10, 0, 1, 3, 5, 10, 0, 1, -3, 6, 10], dtype=float
    )
    grads = rmg.calc_grad_across_cell_corners(values)

    assert_array_equal(
        grads,
        (1.0 / np.sqrt(2.0))
        * np.array(
            [
                [2.0, -1.0, -1.0, 2.0],
                [2.0, -2.0, -2.0, 3.0],
                [4.0, -3.0, -3.0, 4.0],
                [-4.0, -1.0, -1.0, 2.0],
                [3.0, -2.0, -2.0, 3.0],
                [5.0, -8.0, -2.0, 5.0],
            ]
        ),
    )


def test_with_out_keyword():
    """Test with out keyword."""
    rmg = RasterModelGrid((4, 5))
    values_at_nodes = np.arange(20.0)
    out = np.empty((1, 4))
    rtn = rmg.calc_grad_across_cell_corners(values_at_nodes, 5, out=out)
    assert rtn is out
    assert_array_equal(out, np.array([[6.0, 4.0, -6.0, -4.0]]) / np.sqrt(2))

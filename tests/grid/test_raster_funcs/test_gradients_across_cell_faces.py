import numpy as np
from numpy.testing import assert_array_equal

from landlab import RasterModelGrid
from landlab.grid.raster_gradients import calc_grad_across_cell_faces


def test_scalar_arg():
    """Test with a scalar arg for faces."""
    grid = RasterModelGrid((4, 5))
    values_at_nodes = np.arange(20.0)
    grads = calc_grad_across_cell_faces(grid, values_at_nodes, 0)
    assert_array_equal(grads, np.array([[1.0, 5.0, -1.0, -5.0]]))


def test_iterable():
    """Test with an iterable arg for faces."""
    grid = RasterModelGrid((4, 5))
    values_at_nodes = np.arange(20.0)
    grads = grid.calc_grad_across_cell_faces(values_at_nodes, [0, 4])
    assert_array_equal(
        grads, np.array([[1.0, 5.0, -1.0, -5.0], [1.0, 5.0, -1.0, -5.0]])
    )


def test_with_no_cell_id_arg():
    """Test without an arg for faces."""
    grid = RasterModelGrid((4, 5))
    values = np.array(
        [0, 1, 3, 6, 10, 0, 1, 3, 6, 10, 0, 1, 3, 5, 10, 0, 1, -3, 6, 10], dtype=float
    )
    grads = grid.calc_grad_across_cell_faces(values)

    assert_array_equal(
        grads,
        np.array(
            [
                [2.0, 0.0, -1.0, 0.0],
                [3.0, 0.0, -2.0, 0.0],
                [4.0, -1.0, -3.0, 0.0],
                [2.0, 0.0, -1.0, 0.0],
                [2.0, -6.0, -2.0, 0.0],
                [5.0, 1.0, -2.0, 1.0],
            ]
        ),
    )


def test_with_out_keyword():
    """Test using the out keyword."""
    grid = RasterModelGrid((4, 5))
    values_at_nodes = np.arange(20.0)
    out = np.empty((1, 4))
    rtn = grid.calc_grad_across_cell_faces(values_at_nodes, 5, out=out)
    assert rtn is out
    assert_array_equal(out, np.array([[1.0, 5.0, -1.0, -5.0]]))

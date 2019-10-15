import numpy as np
from numpy.testing import assert_array_equal

from landlab import RasterModelGrid


def test_non_unit_spacing():
    """Test on a grid with non-unit spacing."""
    grid = RasterModelGrid((4, 5), xy_spacing=(2, 5))
    values_at_nodes = np.arange(20.0)
    grads = grid.calc_grad_at_link(values_at_nodes)
    assert_array_equal(
        grads,
        np.array(
            [
                0.5,
                0.5,
                0.5,
                0.5,
                1.0,
                1.0,
                1.0,
                1.0,
                1.0,
                0.5,
                0.5,
                0.5,
                0.5,
                1.0,
                1.0,
                1.0,
                1.0,
                1.0,
                0.5,
                0.5,
                0.5,
                0.5,
                1.0,
                1.0,
                1.0,
                1.0,
                1.0,
                0.5,
                0.5,
                0.5,
                0.5,
            ],
            dtype=float,
        ),
    )
    diffs = grid.calc_diff_at_link(values_at_nodes)
    assert_array_equal(
        diffs,
        np.array(
            [
                1,
                1,
                1,
                1,
                5,
                5,
                5,
                5,
                5,
                1,
                1,
                1,
                1,
                5,
                5,
                5,
                5,
                5,
                1,
                1,
                1,
                1,
                5,
                5,
                5,
                5,
                5,
                1,
                1,
                1,
                1,
            ],
            dtype=float,
        ),
    )


def test_out_array():
    """Test using the out keyword."""
    grid = RasterModelGrid((4, 5))
    values_at_nodes = np.arange(20.0)
    grads = np.empty(31)
    rtn_grads = grid.calc_grad_at_link(values_at_nodes, out=grads)
    assert_array_equal(
        grads,
        np.array(
            [
                1,
                1,
                1,
                1,
                5,
                5,
                5,
                5,
                5,
                1,
                1,
                1,
                1,
                5,
                5,
                5,
                5,
                5,
                1,
                1,
                1,
                1,
                5,
                5,
                5,
                5,
                5,
                1,
                1,
                1,
                1,
            ],
            dtype=float,
        ),
    )
    assert rtn_grads is grads

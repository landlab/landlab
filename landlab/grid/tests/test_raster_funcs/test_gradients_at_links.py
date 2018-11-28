import numpy as np
from numpy.testing import assert_array_equal

from landlab import RasterModelGrid


def test_unit_spacing():
    """Test on a grid with unit spacing."""
    grid = RasterModelGrid((4, 5))
    values_at_nodes = np.arange(20.)
    grads = grid.calc_grad_at_link(values_at_nodes)
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
    diffs = grid.calculate_diff_at_links(values_at_nodes)
    assert_array_equal(grads, diffs)


def test_non_unit_spacing():
    """Test on a grid with non-unit spacing."""
    grid = RasterModelGrid((4, 5), spacing=(5, 2))
    values_at_nodes = np.arange(20.)
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
                1.,
                1.,
                1.,
                0.5,
                0.5,
                0.5,
                0.5,
                1.,
                1.,
                1.,
                1.,
                1.,
                0.5,
                0.5,
                0.5,
                0.5,
                1.,
                1.,
                1.,
                1.,
                1.,
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
    values_at_nodes = np.arange(20.)
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


def test_diff_out_array():
    """Test that return array is the out keyword."""
    grid = RasterModelGrid((4, 5))
    values_at_nodes = np.arange(20.)
    diff = np.empty(31)
    rtn_diff = grid.calculate_diff_at_links(values_at_nodes, out=diff)
    assert_array_equal(
        diff,
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
    assert rtn_diff is diff

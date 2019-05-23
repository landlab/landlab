import numpy as np
from numpy.testing import assert_array_equal

from landlab import RasterModelGrid

_GRIDS = {}


def setup_grids():
    """Set up test grids with unit and non-unit spacing."""
    _GRIDS.update(
        {
            "unit": RasterModelGrid((4, 5)),
            "non_unit": RasterModelGrid((4, 5), xy_spacing=2.0),
            "non_square": RasterModelGrid((4, 5), xy_spacing=(5, 2)),
        }
    )


def test_unit_spacing():
    """Test with a grid with unit spacing."""
    rmg = RasterModelGrid((4, 5))
    values_at_nodes = np.arange(20)
    grads = rmg.calc_grad_at_link(values_at_nodes)[rmg.active_links]

    assert_array_equal(
        grads,
        np.array(
            [
                5.0,
                5.0,
                5.0,
                1.0,
                1.0,
                1.0,
                1.0,
                5.0,
                5.0,
                5.0,
                1.0,
                1.0,
                1.0,
                1.0,
                5.0,
                5.0,
                5.0,
            ]
        ),
    )

    diffs = rmg.calc_diff_at_link(values_at_nodes)[rmg.active_links]
    assert_array_equal(grads, diffs)


def test_non_unit_spacing():
    """Test with a grid with non-unit spacing."""
    rmg = RasterModelGrid((4, 5), xy_spacing=(2, 5))
    values_at_nodes = np.arange(20)

    grads = rmg.calc_grad_at_link(values_at_nodes)[rmg.active_links]
    assert_array_equal(
        grads,
        np.array(
            [
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
                0.5,
                0.5,
                0.5,
                0.5,
                1.0,
                1.0,
                1.0,
            ]
        ),
    )
    diffs = rmg.calc_diff_at_link(values_at_nodes)[rmg.active_links]
    assert_array_equal(
        diffs,
        np.array(
            [
                5.0,
                5.0,
                5.0,
                1.0,
                1.0,
                1.0,
                1.0,
                5.0,
                5.0,
                5.0,
                1.0,
                1.0,
                1.0,
                1.0,
                5.0,
                5.0,
                5.0,
            ]
        ),
    )


def test_out_array():
    """Test using the out keyword."""
    rmg = RasterModelGrid((4, 5), xy_spacing=(2, 5))
    values_at_nodes = np.arange(20)

    output_array = np.empty(rmg.number_of_links)
    rtn_array = rmg.calc_grad_at_link(values_at_nodes, out=output_array)
    assert_array_equal(
        rtn_array[rmg.active_links],
        np.array(
            [
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
                0.5,
                0.5,
                0.5,
                0.5,
                1.0,
                1.0,
                1.0,
            ]
        ),
    )
    assert rtn_array is output_array


def test_diff_out_array():
    """Test returned array is the same as that passed as out keyword."""
    rmg = RasterModelGrid((4, 5))
    values = np.arange(20)
    diff = np.empty(rmg.number_of_links)
    rtn_diff = rmg.calc_diff_at_link(values, out=diff)
    assert_array_equal(
        diff[rmg.active_links],
        np.array([5, 5, 5, 1, 1, 1, 1, 5, 5, 5, 1, 1, 1, 1, 5, 5, 5]),
    )
    assert rtn_diff is diff

import numpy as np
from numpy.testing import assert_array_equal
from nose import with_setup
try:
    from nose.tools import assert_is
except ImportError:
    from landlab.testing.tools import assert_is

from landlab import RasterModelGrid


_GRIDS = {}


def setup_grids():
    """Set up test grids with unit and non-unit spacing."""
    _GRIDS.update({
        'unit': RasterModelGrid(4, 5),
        'non_unit': RasterModelGrid(4, 5, 2.),
        'non_square': RasterModelGrid((4, 5), spacing=(5, 2)),
    })


@with_setup(setup_grids)
def test_unit_spacing():
    """Test with a grid with unit spacing."""
    rmg, values_at_nodes = _GRIDS['unit'], np.arange(20)
    grads = rmg.calc_grad_at_link(values_at_nodes)[rmg.active_links]

    assert_array_equal(grads,
                       np.array([5.0, 5.0, 5.0,
                                 1.0, 1.0, 1.0, 1.0,
                                 5.0, 5.0, 5.0,
                                 1.0, 1.0, 1.0, 1.0,
                                 5.0, 5.0, 5.0,]))

    diffs = rmg.calculate_diff_at_active_links(values_at_nodes)
    assert_array_equal(grads, diffs)


@with_setup(setup_grids)
def test_non_unit_spacing():
    """Test with a grid with non-unit spacing."""
    rmg, values_at_nodes = _GRIDS['non_square'], np.arange(20)

    grads = rmg.calc_grad_at_link(values_at_nodes)[rmg.active_links]
    assert_array_equal(grads,
                       np.array([1.0, 1.0, 1.0,
                                 0.5, 0.5, 0.5, 0.5,
                                 1.0, 1.0, 1.0,
                                 0.5, 0.5, 0.5, 0.5,
                                 1.0, 1.0, 1.0]))
    diffs = rmg.calculate_diff_at_active_links(values_at_nodes)
    assert_array_equal(diffs,
                       np.array([5.0, 5.0, 5.0,
                                 1.0, 1.0, 1.0, 1.0,
                                 5.0, 5.0, 5.0,
                                 1.0, 1.0, 1.0, 1.0,
                                 5.0, 5.0, 5.0,]))


@with_setup(setup_grids)
def test_out_array():
    """Test using the out keyword."""
    rmg, values_at_nodes = _GRIDS['non_square'], np.arange(20)

    output_array = np.empty(rmg.number_of_links)
    rtn_array = rmg.calc_grad_at_link(values_at_nodes, out=output_array)
    assert_array_equal(rtn_array[rmg.active_links],
                       np.array([1.0, 1.0, 1.0,
                                 0.5, 0.5, 0.5, 0.5,
                                 1.0, 1.0, 1.0,
                                 0.5, 0.5, 0.5, 0.5,
                                 1.0, 1.0, 1.0]))
    assert_is(rtn_array, output_array)


@with_setup(setup_grids)
def test_diff_out_array():
    """Test returned array is the same as that passed as out keyword."""
    rmg = RasterModelGrid(4, 5)
    values = np.arange(20)
    diff = np.empty(17)
    rtn_diff = rmg.calculate_diff_at_active_links(values, out=diff)
    assert_array_equal(
        diff,
        np.array([5, 5, 5,
                  1, 1, 1, 1,
                  5, 5, 5,
                  1, 1, 1, 1,
                  5, 5, 5]))
    assert_is(rtn_diff, diff)

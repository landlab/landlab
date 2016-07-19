import numpy as np
from numpy.testing import assert_array_equal
from nose import with_setup
from nose.tools import assert_equal
try:
    from nose.tools import assert_is
except ImportError:
    from landlab.testing.tools import assert_is


def setup_unit_grid():
    from landlab import RasterModelGrid
    globals().update({
        'rmg': RasterModelGrid(4, 5),
        'values_at_nodes':  np.arange(20.),
    })


def setup_3x3_grid():
    from landlab import RasterModelGrid
    globals().update({
        'rmg_3x3': RasterModelGrid(3, 3),
    })


def setup_grid():
    from landlab import RasterModelGrid
    globals().update({
        'rmg': RasterModelGrid(4, 5, 2.),
        'values_at_nodes':  np.arange(20.),
    })


@with_setup(setup_unit_grid)
def test_scalar_arg():
    grad = rmg._calc_steepest_descent_across_cell_corners(
        values_at_nodes, 0)
    assert_equal(grad, - 6. / np.sqrt(2.))


@with_setup(setup_unit_grid)
def test_iterable():
    grads = rmg._calc_steepest_descent_across_cell_corners(values_at_nodes,
                                                               [0, 4])
    assert_array_equal(grads, - np.array([6., 6.]) / np.sqrt(2))


@with_setup(setup_unit_grid)
def test_scalar_arg_with_faces_ids():
    values = np.array([0, 1,  3, 6, 10,
                       0, 1,  3, 6, 10,
                       0, 1,  3, 5, 10,
                       0, 1, -3, 6, 10], dtype=float)
    (grad, node) = rmg._calc_steepest_descent_across_cell_corners(
        values, (0, 4), return_node=True)
    assert_array_equal(grad, - np.array([1, 2]) / np.sqrt(2.))
    assert_array_equal(node, [10, 16])


@with_setup(setup_3x3_grid)
def test_node_in_direction_of_max():
    for diagonal_id in [0, 2, 6, 8]:
        values = np.zeros(9)
        values[diagonal_id] = -1
        (_, node) = rmg_3x3._calc_steepest_descent_across_cell_corners(
            values, 0, return_node=True)
        assert_array_equal(node, diagonal_id)


@with_setup(setup_3x3_grid)
def test_node_in_direction_of_max_with_ties():
    values = np.zeros(9)
    (_, node) = rmg_3x3._calc_steepest_descent_across_cell_corners(
        values, 0, return_node=True)
    assert_array_equal(node, 8)

    values = np.zeros(9)
    values[8] = 1
    (_, node) = rmg_3x3._calc_steepest_descent_across_cell_corners(
        values, 0, return_node=True)
    assert_array_equal(node, 6)

    values = np.zeros(9)
    values[[8, 6]] = 1
    (_, node) = rmg_3x3._calc_steepest_descent_across_cell_corners(
        values, 0, return_node=True)
    assert_array_equal(node, 0)

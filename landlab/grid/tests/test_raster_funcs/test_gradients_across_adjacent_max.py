import numpy as np
from numpy.testing import assert_array_equal
from nose import with_setup
try:
    from nose.tools import assert_is
except ImportError:
    from landlab.testing.tools import assert_is
from nose.tools import assert_equal

from landlab.grid.raster_steepest_descent import (
    _calc_steepest_descent_across_adjacent_cells)

def setup_unit_grid():
    """Set up a test grid with unit spacing."""
    from landlab import RasterModelGrid
    globals()['rmg'] = RasterModelGrid(4, 5)
    globals()['values_at_nodes'] = np.arange(20, dtype=float)


def setup_non_unit_grid():
    """Set up a test grid with non-unit spacing."""
    from landlab import RasterModelGrid
    globals()['rmg'] = RasterModelGrid(4, 5, 2)
    globals()['values_at_nodes'] = np.arange(20, dtype=float)


def setup_3x3_grid():
    """Set up a grid of 3 rows and 3 columns."""
    from landlab import RasterModelGrid
    globals().update({
        'rmg_3x3': RasterModelGrid(3, 3),
        'values_at_nodes': np.flipud(np.array([6, 7, 8,
                                               3, 4, 5,
                                               0, 1, 2], dtype=float))
    })


@with_setup(setup_3x3_grid)
def test_scalar_arg():
    """Test scalar arg for nodes."""
    grad = _calc_steepest_descent_across_adjacent_cells(
        rmg_3x3, values_at_nodes, 0)
    assert_equal(grad, -3.)

    grad = _calc_steepest_descent_across_adjacent_cells(
        rmg_3x3, values_at_nodes, 0, method='d8')
    assert_equal(grad, -3.)

    values_at_nodes[2] = -10
    grad = _calc_steepest_descent_across_adjacent_cells(
        rmg_3x3, values_at_nodes, 0, method='d8')
    assert_equal(grad, - (4 + 10) / np.sqrt(2.))


@with_setup(setup_unit_grid)
def test_iterable():
    """Test iterable arg for nodes."""
    grad = _calc_steepest_descent_across_adjacent_cells(
        rmg, values_at_nodes, [0, 4])
    assert_array_equal(grad, [-5., -5.])


@with_setup(setup_unit_grid)
def test_scalar_arg_with_links():
    values = np.array([0, 1,  3, 6, 10,
                       0, 1,  3, 6, 10,
                       0, 1,  3, 5, 10,
                       0, 1, -3, 6, 10, ], dtype=float)
    (grad, node) = _calc_steepest_descent_across_adjacent_cells(
        rmg, values, (0, 4), return_node=True)
    assert_array_equal(grad, [-1, -6])
    assert_array_equal(node, [5, 17])

    values_at_nodes[2] = -10
    (grad, node) = _calc_steepest_descent_across_adjacent_cells(
        rmg, values_at_nodes, 0, method='d8', return_node=True)
    assert_equal(grad, - (6 + 10) / np.sqrt(2.))
    assert_equal(node, 2)


@with_setup(setup_unit_grid)
def test_node_id_in_direction_of_max():
    values = np.array([-1, 1,  3, 6, 10,
                       0, 1,  3, 6, 10,
                       0, 1,  3, 5, 10,
                       0, 1, -3, 6, 10, ], dtype=float)
    (_, node_ids) = _calc_steepest_descent_across_adjacent_cells(
        rmg, values, (0, 4), return_node=True)
    assert_array_equal(node_ids, [5, 17])

    (grads,
     node_ids) = _calc_steepest_descent_across_adjacent_cells(
        rmg, values, (0, 4), method='d8', return_node=True)
    assert_array_equal(node_ids, [0, 17])


@with_setup(setup_3x3_grid)
def test_node_in_direction_of_max():
    for node_id in [0, 1, 2, 3, 5, 6, 7, 8]:
        values = np.zeros(9)
        values[node_id] = -1
        (_, node) = _calc_steepest_descent_across_adjacent_cells(
            rmg_3x3, values, 0, return_node=True, method='d8')
        assert_array_equal(node, node_id)


@with_setup(setup_3x3_grid)
def test_node_in_direction_of_max_with_ties():
    values = np.zeros(9)
    (_, node) = _calc_steepest_descent_across_adjacent_cells(
        rmg_3x3, values, 0, return_node=True, method='d8')
    assert_array_equal(node, 5)

    for (node_id, expected) in zip([5, 7, 3, 1, 8, 6, 0],
                                   [7, 3, 1, 8, 6, 0, 2]):
        values[node_id] = 1
        (_, node) = _calc_steepest_descent_across_adjacent_cells(
            rmg_3x3, values, 0, return_node=True, method='d8')
        assert_array_equal(node, expected)

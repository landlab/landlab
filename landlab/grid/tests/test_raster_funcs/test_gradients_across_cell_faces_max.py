import numpy as np
from numpy.testing import assert_array_equal
from nose import with_setup
from nose.tools import assert_equal
try:
    from nose.tools import assert_is
except ImportError:
    from landlab.testing.tools import assert_is

from landlab.grid.raster_steepest_descent import (
    _calc_steepest_descent_across_cell_faces)


_GRID = None
_VALUES_AT_NODES = None


def setup_unit_grid():
    from landlab import RasterModelGrid
    global _GRID, _VALUES_AT_NODES
    _GRID = RasterModelGrid(4, 5)
    _VALUES_AT_NODES = np.arange(20.)


def setup_3x3_grid():
    from landlab import RasterModelGrid
    global rmg_3x3
    rmg_3x3 = RasterModelGrid(3, 3)


@with_setup(setup_unit_grid)
def test_scalar_arg():
    grad = _calc_steepest_descent_across_cell_faces(
        _GRID, _VALUES_AT_NODES, 0)
    assert_equal(grad, -5.)

    grad = _GRID._calc_steepest_descent_across_cell_faces(
        _VALUES_AT_NODES, 0)
    assert_equal(grad, -5.)


@with_setup(setup_unit_grid)
def test_iterable():
    grad = _calc_steepest_descent_across_cell_faces(
        _GRID, _VALUES_AT_NODES, [0, 4])
    assert_array_equal(grad, [-5., -5.])

    grad = _GRID._calc_steepest_descent_across_cell_faces(
        _VALUES_AT_NODES, [0, 4])
    assert_array_equal(grad, [-5., -5.])


@with_setup(setup_unit_grid)
def test_scalar_arg_with_return_node():
    values = np.array([0, 1,  3, 6, 10,
                       0, 1,  3, 6, 10,
                       0, 1,  3, 5, 10,
                       0, 1, -3, 6, 10, ], dtype=float)
    (grad, node) = _GRID._calc_steepest_descent_across_cell_faces(
        values, (0, 4), return_node=True)
    assert_array_equal(grad, [-1, -6])
    assert_array_equal(node, [5, 17])


@with_setup(setup_3x3_grid)
def test_node_in_direction_of_max():
    for neighbor_id in [1, 3, 5, 7]:
        values = np.zeros(9)
        values[neighbor_id] = -1
        (_, node) = rmg_3x3._calc_steepest_descent_across_cell_faces(
            values, 0, return_node=True)
        assert_array_equal(node, neighbor_id)


@with_setup(setup_3x3_grid)
def test_node_in_direction_of_max_with_ties():
    values = np.zeros(9)
    (_, node) = rmg_3x3._calc_steepest_descent_across_cell_faces(
        values, 0, return_node=True)
    assert_array_equal(node, 5)

    values = np.zeros(9)
    values[5] = 1
    (_, node) = rmg_3x3._calc_steepest_descent_across_cell_faces(
        values, 0, return_node=True)
    assert_array_equal(node, 7)

    values = np.zeros(9)
    values[[5, 7]] = 1
    (_, node) = rmg_3x3._calc_steepest_descent_across_cell_faces(
        values, 0, return_node=True)
    assert_array_equal(node, 3)

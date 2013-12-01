import numpy as np
from numpy.testing import assert_array_equal
from nose import with_setup
from nose.tools import assert_is, assert_equal

from landlab.grid import raster_funcs as rfuncs


def setup_unit_grid():
    from landlab import RasterModelGrid
    globals()['rmg'] = RasterModelGrid(4, 5)
    globals()['values_at_nodes'] = np.arange(20)


def setup_non_unit_grid():
    from landlab import RasterModelGrid
    globals()['rmg'] = RasterModelGrid(4, 5, 2)
    globals()['values_at_nodes'] = np.arange(20)


@with_setup(setup_unit_grid)
def test_scalar_arg():
    grad = rfuncs.calculate_max_gradient_across_adjacent_cells(
        rmg, values_at_nodes, 0)
    assert_equal(grad, 5.)

    grad = rfuncs.calculate_max_gradient_across_adjacent_cells(
        rmg, values_at_nodes, 0, method='d8')
    assert_equal(grad, 5.)

    values_at_nodes[2] = -10
    grad = rfuncs.calculate_max_gradient_across_adjacent_cells(
        rmg, values_at_nodes, 0, method='d8')
    assert_equal(grad, 11)


@with_setup(setup_unit_grid)
def test_iterable():
    grad = rfuncs.calculate_max_gradient_across_adjacent_cells(
        rmg, values_at_nodes, [0, 4])
    assert_array_equal(grad, [5., 5.])


@with_setup(setup_unit_grid)
def test_scalar_arg_with_links():
    values = np.array([0, 1,  3, 6, 10,
                       0, 1,  3, 6, 10,
                       0, 1,  3, 5, 10,
                       0, 1, -3, 6, 10,])
    (grad, face) = rfuncs.calculate_max_gradient_across_adjacent_cells(
        rmg, values, (0, 4), return_face=True)
    assert_array_equal(grad, [1, 6])
    assert_array_equal(face, [1, 2])

    link_ids = rfuncs.active_link_id_of_cell_neighbor(rmg, face, [0, 4])
    assert_array_equal(link_ids, [9, 7])

    node_ids = rfuncs.node_id_of_cell_neighbor(rmg, face, (0, 4))
    assert_array_equal(node_ids, [5, 17])

    values_at_nodes[2] = -10
    (grad, face) = rfuncs.calculate_max_gradient_across_adjacent_cells(
        rmg, values_at_nodes, 0, method='d8', return_face=True)
    assert_equal(grad, 11)
    assert_equal(face, 4)



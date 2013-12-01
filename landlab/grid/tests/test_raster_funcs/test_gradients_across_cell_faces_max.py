import numpy as np
from numpy.testing import assert_array_equal
from nose import with_setup
from nose.tools import assert_is, assert_equal

from landlab.grid import raster_funcs as rfuncs


_GRID = None
_VALUES_AT_NODES = None


def setup_unit_grid():
    from landlab import RasterModelGrid
    global _GRID, _VALUES_AT_NODES
    _GRID = RasterModelGrid(4, 5)
    _VALUES_AT_NODES = np.arange(20)


@with_setup(setup_unit_grid)
def test_scalar_arg():
    grad = rfuncs.calculate_max_gradient_across_cell_faces(_GRID,
                                                           _VALUES_AT_NODES, 0)
    assert_equal(grad, 5.)

    grad = _GRID.calculate_max_gradient_across_cell_faces(_VALUES_AT_NODES, 0)
    assert_equal(grad, 5.)


@with_setup(setup_unit_grid)
def test_iterable():
    grad = rfuncs.calculate_max_gradient_across_cell_faces(
        _GRID, _VALUES_AT_NODES, [0, 4])
    assert_array_equal(grad, [5., 5.])

    grad = _GRID.calculate_max_gradient_across_cell_faces(
        _VALUES_AT_NODES, [0, 4])
    assert_array_equal(grad, [5., 5.])


@with_setup(setup_unit_grid)
def test_scalar_arg_with_links():
    values = np.array([0, 1,  3, 6, 10,
                       0, 1,  3, 6, 10,
                       0, 1,  3, 5, 10,
                       0, 1, -3, 6, 10,])
    (grad, face) = _GRID.calculate_max_gradient_across_cell_faces(
        values, (0, 4), return_face=True)
    assert_array_equal(grad, [1, 6])
    assert_array_equal(face, [1, 2])

    link_ids = rfuncs.active_link_id_of_cell_neighbor(_GRID, face, [0, 4])
    assert_array_equal(link_ids, [9, 7])

    node_ids = rfuncs.node_id_of_cell_neighbor(_GRID, face, (0, 4))
    assert_array_equal(node_ids, [5, 17])

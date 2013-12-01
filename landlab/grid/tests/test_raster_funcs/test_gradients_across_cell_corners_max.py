import numpy as np
from numpy.testing import assert_array_equal
from nose import with_setup
from nose.tools import assert_is, assert_equal

from landlab.grid import raster_funcs as rfuncs


def setup_unit_grid():
    from landlab import RasterModelGrid
    globals().update({
        'rmg': RasterModelGrid(4, 5),
        'values_at_nodes':  np.arange(20.),
    })


def setup_grid():
    from landlab import RasterModelGrid
    globals().update({
        'rmg': RasterModelGrid(4, 5, 2.),
        'values_at_nodes':  np.arange(20.),
    })


@with_setup(setup_unit_grid)
def test_scalar_arg():
    grad = rmg.calculate_max_gradient_across_cell_corners(values_at_nodes, 0)
    assert_equal(grad, 6. / np.sqrt(2.))


@with_setup(setup_unit_grid)
def test_iterable():
    grads = rmg.calculate_max_gradient_across_cell_corners(values_at_nodes,
                                                           [0, 4])
    assert_array_equal(grads, np.array([6., 6.]) / np.sqrt(2))


@with_setup(setup_unit_grid)
def test_scalar_arg_with_faces_ids():
    values = np.array([0, 1,  3, 6, 10,
                       0, 1,  3, 6, 10,
                       0, 1,  3, 5, 10,
                       0, 1, -3, 6, 10], dtype=float)
    (grad, face) = rmg.calculate_max_gradient_across_cell_corners(
        values, (0, 4), return_face=True)
    assert_array_equal(grad, np.array([1, 2]) / np.sqrt(2.))
    assert_array_equal(face, [2, 2])

    node_ids = rfuncs.node_id_of_cell_corner(rmg, face, [0, 4])
    assert_array_equal(node_ids, [10, 16])

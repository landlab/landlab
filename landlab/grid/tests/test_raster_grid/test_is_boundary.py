import numpy as np
from numpy.testing import assert_array_equal
from nose.tools import with_setup, assert_true

from landlab import FIXED_GRADIENT_BOUNDARY
from landlab import RasterModelGrid


def setup_grid():
    globals().update({
        'rmg': RasterModelGrid(4, 5)
    })


@with_setup(setup_grid)
def test_id_as_int():
    assert_true(rmg.node_is_boundary(0))


@with_setup(setup_grid)
def test_id_as_small_list():
    assert_array_equal(rmg.node_is_boundary([0]), np.array([True]))


@with_setup(setup_grid)
def test_id_as_array():
    assert_array_equal(
        rmg.node_is_boundary(np.arange(20)),
        np.array([True,  True,  True,  True,  True,
                  True, False, False, False,  True,
                  True, False, False, False,  True,
                  True,  True,  True,  True,  True], dtype=bool))


@with_setup(setup_grid)
def test_id_as_list():
    assert_array_equal(rmg.node_is_boundary([8, 9]), np.array([False, True]))


@with_setup(setup_grid)
def test_boundary_flag():
    rmg.status_at_node[0] = FIXED_GRADIENT_BOUNDARY
    assert_array_equal(
        rmg.node_is_boundary(np.arange(20)),
        np.array([True,  True,  True,  True,  True,
                  True, False, False, False,  True,
                  True, False, False, False,  True,
                  True,  True,  True,  True,  True], dtype=bool))

    assert_array_equal(
        rmg.node_is_boundary(np.arange(20), boundary_flag=FIXED_GRADIENT_BOUNDARY),
        np.array([True, False, False, False, False,
                  False, False, False, False, False,
                  False, False, False, False, False,
                  False, False, False, False, False], dtype=bool))

import numpy as np
from numpy.testing import assert_array_equal
from nose.tools import with_setup, assert_true, assert_equal, assert_not_equal

from landlab import FIXED_GRADIENT_BOUNDARY, CLOSED_BOUNDARY, INACTIVE_LINK, \
    FIXED_LINK
from landlab import RasterModelGrid


def setup_grid():
    globals().update({
        'rmg': RasterModelGrid((4, 5))
    })


@with_setup(setup_grid)
def test_link_update_with_nodes_closed():
    rmg.status_at_node[rmg.nodes_at_bottom_edge] = CLOSED_BOUNDARY
    inactive_array = np.array([INACTIVE_LINK, ] * 5)
    assert_array_equal(rmg.status_at_link[4:9], inactive_array)


@with_setup(setup_grid)
def test_link_update_with_nodes_fixed_grad():
    rmg.status_at_node[rmg.nodes_at_bottom_edge] = FIXED_GRADIENT_BOUNDARY
    fixed_array = np.array([FIXED_LINK, ] * 3)
    assert_array_equal(rmg.status_at_link[5:8], fixed_array)


@with_setup(setup_grid)
def test_bc_set_code_init():
    assert_equal(rmg.bc_set_code, 0)


@with_setup(setup_grid)
def test_bc_set_code_change():
    rmg.status_at_node[rmg.nodes_at_bottom_edge] = CLOSED_BOUNDARY
    assert_not_equal(rmg.bc_set_code, 0)

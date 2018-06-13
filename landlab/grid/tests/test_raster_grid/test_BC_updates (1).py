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


def test_issue_428_a():
    """Issue #428"""
    grid = RasterModelGrid((4, 4))
    grid.set_closed_boundaries_at_grid_edges(True, True, True, True)

    assert_equal(grid.status_at_node[1], 4)
    assert_equal(grid.status_at_link[4], 4)
    assert_array_equal(grid.active_link_dirs_at_node[1],[0, 0, 0, 0])

    grid.status_at_node[1] = 1
    assert_equal(grid.status_at_link[4], 0)
    assert_array_equal(grid.active_link_dirs_at_node[1], [0, -1, 0, 0])


def test_issue_428_b():
    """Issue #428"""
    grid = RasterModelGrid((4, 4))

    z = np.ones(grid.number_of_nodes)
    z[grid.nodes_at_bottom_edge] = -9999.
    z[grid.nodes_at_left_edge] = -9999.
    z[grid.nodes_at_top_edge] = -9999.
    z[grid.nodes_at_right_edge] = -9999.
    z[1] = .5

    assert_array_equal(grid.active_link_dirs_at_node[1], [0, -1, 0, 0])

    grid.set_watershed_boundary_condition(z)
    assert_array_equal(grid.active_link_dirs_at_node[1], [0, -1, 0, 0])


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

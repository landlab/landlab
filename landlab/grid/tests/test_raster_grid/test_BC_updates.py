import numpy as np
from numpy.testing import assert_array_equal

from landlab import (
    CLOSED_BOUNDARY,
    FIXED_GRADIENT_BOUNDARY,
    FIXED_LINK,
    INACTIVE_LINK,
    RasterModelGrid,
)


def test_issue_428_a():
    """Issue #428"""
    grid = RasterModelGrid((4, 4))
    grid.set_closed_boundaries_at_grid_edges(True, True, True, True)

    assert grid.status_at_node[1] == 4
    assert grid.status_at_link[4] == 4
    assert_array_equal(grid.active_link_dirs_at_node[1], [0, 0, 0, 0])

    grid.status_at_node[1] = 1
    assert grid.status_at_link[4] == 0
    assert_array_equal(grid.active_link_dirs_at_node[1], [0, -1, 0, 0])


def test_issue_428_b():
    """Issue #428"""
    grid = RasterModelGrid((4, 4))

    z = np.ones(grid.number_of_nodes)
    z[grid.nodes_at_bottom_edge] = -9999.0
    z[grid.nodes_at_left_edge] = -9999.0
    z[grid.nodes_at_top_edge] = -9999.0
    z[grid.nodes_at_right_edge] = -9999.0
    z[1] = 0.5

    assert_array_equal(grid.active_link_dirs_at_node[1], [0, -1, 0, 0])

    grid.set_watershed_boundary_condition(z)
    assert_array_equal(grid.active_link_dirs_at_node[1], [0, -1, 0, 0])


def test_link_update_with_nodes_closed():
    rmg = RasterModelGrid((4, 5))
    rmg.status_at_node[rmg.nodes_at_bottom_edge] = CLOSED_BOUNDARY
    inactive_array = np.array([INACTIVE_LINK] * 5)
    assert_array_equal(rmg.status_at_link[4:9], inactive_array)


def test_link_update_with_nodes_fixed_grad():
    rmg = RasterModelGrid((4, 5))
    rmg.status_at_node[rmg.nodes_at_bottom_edge] = FIXED_GRADIENT_BOUNDARY
    fixed_array = np.array([FIXED_LINK] * 3)
    assert_array_equal(rmg.status_at_link[5:8], fixed_array)


def test_bc_set_code_init():
    rmg = RasterModelGrid((4, 5))
    assert rmg.bc_set_code == 0


def test_bc_set_code_change():
    rmg = RasterModelGrid((4, 5))
    rmg.status_at_node[rmg.nodes_at_bottom_edge] = CLOSED_BOUNDARY
    assert rmg.bc_set_code != 0

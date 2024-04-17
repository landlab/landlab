"""
Created on Wed Jan 13 15:34:08 2021

@author: laure
"""

import numpy as np

from landlab import HexModelGrid
from landlab import RadialModelGrid
from landlab import RasterModelGrid
from landlab import VoronoiDelaunayGrid
from landlab.utils.window_statistic import calculate_window_statistic


def test_hex_grid():
    grid = HexModelGrid((5, 4), spacing=10)
    grid.status_at_node[grid.status_at_node == 1] = grid.BC_NODE_IS_CLOSED
    z = grid.add_zeros("topographic__elevation", at="node")
    z += np.arange(len(z))
    relief = calculate_window_statistic(
        grid,
        "topographic__elevation",
        np.ptp,
        search_radius=15,
        calc_on_closed_nodes=False,
    )

    check = [
        np.nan,
        np.nan,
        np.nan,
        np.nan,
        np.nan,
        6.0,
        7.0,
        7.0,
        np.nan,
        np.nan,
        11.0,
        12.0,
        12.0,
        11.0,
        np.nan,
        np.nan,
        7.0,
        7.0,
        6.0,
        np.nan,
        np.nan,
        np.nan,
        np.nan,
        np.nan,
    ]

    np.testing.assert_equal(relief, check)


def test_calc_with_interior_closed_nodes():
    grid = RasterModelGrid((5, 5), xy_spacing=10.0)
    grid.status_at_node[[7, 8, 12, 13]] = grid.BC_NODE_IS_CLOSED
    z = grid.add_zeros("topographic__elevation", at="node")
    z += np.arange(len(z))
    relief = calculate_window_statistic(
        grid,
        "topographic__elevation",
        np.ptp,
        search_radius=15,
        calc_on_closed_nodes=False,
    )

    check = [
        6.0,
        6.0,
        5.0,
        7.0,
        6.0,
        11.0,
        11.0,
        np.nan,
        np.nan,
        11.0,
        11.0,
        12.0,
        np.nan,
        np.nan,
        10.0,
        11.0,
        12.0,
        12.0,
        10.0,
        10.0,
        6.0,
        7.0,
        7.0,
        7.0,
        6.0,
    ]

    np.testing.assert_equal(relief, check)


def test_voronoi_grid():
    x = [0, 0.001, 0.002, 0.003, 1, 1.001, 1.002, 1.003, 2, 2.001, 2.002, 2.003]
    y = [0, 1, 2, 3, 0, 1, 2, 3, 0, 1, 2, 3]
    grid = VoronoiDelaunayGrid(x, y)
    z = grid.add_zeros("topographic__elevation", at="node")
    z += np.arange(len(z))
    relief = calculate_window_statistic(
        grid, "topographic__elevation", np.ptp, search_radius=1.5
    )

    check = [4.0, 5.0, 4.0, 7.0, 8.0, 7.0, 7.0, 8.0, 7.0, 4.0, 5.0, 4.0]

    np.testing.assert_equal(relief, check)


def test_radial_grid():
    grid = RadialModelGrid(2, 6, 1)
    grid.status_at_node[grid.status_at_node == 1] = grid.BC_NODE_IS_CLOSED
    z = grid.add_zeros("topographic__elevation", at="node")
    z += np.arange(len(z))
    relief = calculate_window_statistic(
        grid,
        "topographic__elevation",
        np.ptp,
        search_radius=1.5,
        calc_on_closed_nodes=False,
    )

    check = [
        np.nan,
        np.nan,
        np.nan,
        np.nan,
        np.nan,
        4.0,
        5.0,
        np.nan,
        7.0,
        8.0,
        7.0,
        np.nan,
        5.0,
        4.0,
        np.nan,
        np.nan,
        np.nan,
        np.nan,
        np.nan,
    ]

    np.testing.assert_equal(relief, check)

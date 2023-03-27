import numpy as np
from numpy.testing import assert_array_equal

from landlab import HexModelGrid, RasterModelGrid
from landlab.components.advection import find_upwind_link_at_link


def test_upwind_link_at_link_raster():
    grid = RasterModelGrid((3, 4))
    uwl = find_upwind_link_at_link(grid, 1.0)
    assert_array_equal(
        uwl, [-1, -1, -1, -1, -1, -1, -1, -1, 7, 8, -1, 4, 5, -1, -1, -1, -1]
    )
    uwl = find_upwind_link_at_link(grid, -1.0)
    assert_array_equal(
        uwl, [-1, -1, -1, -1, 11, 12, -1, 8, 9, -1, -1, -1, -1, -1, -1, -1, -1]
    )
    u = np.zeros(grid.number_of_links)
    u[4:6] = -1
    u[7] = -1
    u[8:10] = 1
    u[11:13] = 1
    uwl = find_upwind_link_at_link(grid, u)
    assert_array_equal(
        uwl, [-1, -1, -1, -1, 11, 12, -1, 8, 7, 8, -1, 4, 5, -1, -1, -1, -1]
    )


def test_upwind_link_at_link_hex():

    # Hex horizontal
    grid = HexModelGrid((3, 3))
    uwl = find_upwind_link_at_link(grid, 1.0)
    assert_array_equal(
        uwl, [-1, -1, -1, -1, -1, -1, -1, -1, -1, 8, 9, -1, 4, 3, 6, 5, -1, -1, -1]
    )
    uwl = find_upwind_link_at_link(grid, -1.0)
    assert_array_equal(
        uwl, [-1, -1, -1, 13, 12, 15, 14, -1, 9, 10, -1, -1, -1, -1, -1, -1, -1, -1, -1]
    )
    u = np.zeros(grid.number_of_links)
    u[3:7] = -1
    u[8] = -1
    u[9:11] = 1
    u[12:16] = 1
    uwl = find_upwind_link_at_link(grid, u)
    assert_array_equal(
        uwl, [-1, -1, -1, 13, 12, 15, 14, -1, 9, 8, 9, -1, 4, 3, 6, 5, -1, -1, -1]
    )

    # Hex vertical
    grid = HexModelGrid((3, 3), orientation="vertical")
    uwl = find_upwind_link_at_link(grid, 1.0)
    assert_array_equal(
        uwl, [-1, -1, -1, -1, 7, -1, -1, -1, 3, 2, -1, 14, -1, -1, -1, 10, 9, -1, -1]
    )
    uwl = find_upwind_link_at_link(grid, -1.0)
    assert_array_equal(
        uwl, [-1, -1, 9, 8, -1, -1, -1, 4, -1, 16, 15, -1, -1, -1, 11, -1, -1, -1, -1]
    )
    u = np.zeros(grid.number_of_links)
    u[2:4] = -1
    u[4] = 1
    u[7] = -1
    u[9] = 1
    u[10] = -1
    u[11] = 1
    u[14] = -1
    u[15] = 1
    uwl = find_upwind_link_at_link(grid, u)
    assert_array_equal(
        uwl, [-1, -1, 9, 8, 7, -1, -1, 4, 3, 2, 15, 14, -1, -1, 11, 10, 9, -1, -1]
    )

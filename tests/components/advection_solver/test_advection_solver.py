import numpy as np
import pytest
from numpy.testing import assert_array_almost_equal, assert_array_equal

from landlab import HexModelGrid, RasterModelGrid
from landlab.components.advection import (
    find_upwind_link_at_link,
    upwind_to_local_grad_ratio,
)
from landlab.grid.linkorientation import LinkOrientation


@pytest.mark.parametrize("u", [0.0, 1.0, 10.0, [1.0] * 17])
def test_upwind_link_at_link_raster_positive(u):
    grid = RasterModelGrid((3, 4))
    uwl = find_upwind_link_at_link(grid, u)

    assert_array_equal(
        uwl[grid.vertical_links].reshape((2, 4)),
        [
            [-1, -1, -1, -1],
            [3, 4, 5, 6],
        ],
    )
    assert_array_equal(
        uwl[grid.horizontal_links].reshape((3, 3)),
        [
            [-1, 0, 1],
            [-1, 7, 8],
            [-1, 14, 15],
        ],
    )


@pytest.mark.parametrize("u", [-1.0, -10.0, [-1.0] * 17])
def test_upwind_link_at_link_raster_negative(u):
    grid = RasterModelGrid((3, 4))
    uwl = find_upwind_link_at_link(grid, u)
    assert_array_equal(
        uwl[grid.vertical_links].reshape((2, 4)),
        [
            [10, 11, 12, 13],
            [-1, -1, -1, -1],
        ],
    )
    assert_array_equal(
        uwl[grid.horizontal_links].reshape((3, 3)),
        [
            [1, 2, -1],
            [8, 9, -1],
            [15, 16, -1],
        ],
    )


def test_upwind_link_at_link_raster_mixed():
    grid = RasterModelGrid((3, 4))
    u = np.zeros(grid.number_of_links)
    u[4:6] = -1
    u[7] = -1
    u[8:10] = 1
    u[11:13] = 1
    uwl = find_upwind_link_at_link(grid, u)
    assert_array_equal(
        uwl[grid.vertical_links].reshape((2, 4)),
        [
            [-1, 11, 12, -1],
            [3, 4, 5, 6],
        ],
    )
    assert_array_equal(
        uwl[grid.horizontal_links].reshape((3, 3)),
        [
            [-1, 0, 1],
            [8, 7, 8],
            [-1, 14, 15],
        ],
    )


def test_upwind_link_at_link_hex_horizontal():
    # Hex horizontal
    grid = HexModelGrid((3, 3), orientation="horizontal")
    uwl = find_upwind_link_at_link(grid, 1.0)

    assert_array_equal(
        uwl[grid.orientation_of_link == LinkOrientation.E], [-1, 0, -1, 8, 9, -1, 17]
    )
    assert_array_equal(
        uwl[grid.orientation_of_link == LinkOrientation.NNE],
        [-1, -1, -1, -1, 3, 5],
    )
    assert_array_equal(
        uwl[grid.orientation_of_link == LinkOrientation.NNW], [-1, -1, -1, 4, 6, -1]
    )

    uwl = find_upwind_link_at_link(grid, -1.0)
    assert_array_equal(
        uwl[grid.orientation_of_link == LinkOrientation.E],
        [1, -1, 9, 10, -1, 18, -1],
    )
    assert_array_equal(
        uwl[grid.orientation_of_link == LinkOrientation.NNE],
        [13, 15, -1, -1, -1, -1],
    )
    assert_array_equal(
        uwl[grid.orientation_of_link == LinkOrientation.NNW],
        [-1, 12, 14, -1, -1, -1],
    )

    u = np.zeros(grid.number_of_links)
    u[3:7] = -1
    u[8] = -1
    u[9:11] = 1
    u[12:16] = 1

    expected = np.choose(
        u >= 0, [find_upwind_link_at_link(grid, -1), find_upwind_link_at_link(grid, 1)]
    )
    assert_array_equal(find_upwind_link_at_link(grid, u), expected)


def test_upwind_link_at_link_hex_vertical():
    # Hex vertical
    grid = HexModelGrid((3, 3), orientation="vertical")
    uwl = find_upwind_link_at_link(grid, 1.0)

    assert_array_equal(
        uwl[grid.orientation_of_link == LinkOrientation.ENE], [-1, -1, 3, -1, 10, -1]
    )
    assert_array_equal(
        uwl[grid.orientation_of_link == LinkOrientation.N], [-1, -1, -1, 2, 5, 6, 9]
    )
    assert_array_equal(
        uwl[grid.orientation_of_link == LinkOrientation.ESE], [-1, 7, -1, 14, -1, -1]
    )

    uwl = find_upwind_link_at_link(grid, -1.0)

    assert_array_equal(
        uwl[grid.orientation_of_link == LinkOrientation.ENE], [-1, 8, -1, 15, -1, -1]
    )
    assert_array_equal(
        uwl[grid.orientation_of_link == LinkOrientation.N], [9, 12, 13, 16, -1, -1, -1]
    )
    assert_array_equal(
        uwl[grid.orientation_of_link == LinkOrientation.ESE], [-1, -1, 4, -1, 11, -1]
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

    expected = np.choose(
        u >= 0, [find_upwind_link_at_link(grid, -1), find_upwind_link_at_link(grid, 1)]
    )
    assert_array_equal(find_upwind_link_at_link(grid, u), expected)


def test_upwind_to_local_grad_ratio():
    """
    Predicted upwind_to_local_grad_ratio for u>1

    Link  Local grad  Upwind link  Upwind diff  Ratio
     0 (n/a)
     1 (n/a)
     2 (n/a)
     3 (n/a)
     4    24          none         n/a          1
     5    32          none         n/a          1
     6 (n/a)
     7     9          none         n/a          1
     8    11          7             9           9/11
     9    13          8            11           11/13
    10 (n/a)
    11    56          4            24           24/56
    12    64          5            32           1/2
    13 (n/a)
    14    17          none         n/a          1
    15    19          14           17           17/19
    16    21          15           19           19/21
    17 (n/a)
    18    88          11           56           56/88
    19    96          12           64           64/96
    20 (n/a)
    21 (n/a)
    22 (n/a)
    23 (n/a)
    """
    grid = RasterModelGrid((4, 4))
    v = np.arange(grid.number_of_nodes) ** 2

    upwind_link_at_link = find_upwind_link_at_link(grid, 1.0)
    r = upwind_to_local_grad_ratio(grid, v, upwind_link_at_link)

    diff_at_link = grid.calc_grad_at_link(v)
    expected = diff_at_link[upwind_link_at_link] / diff_at_link

    assert_array_almost_equal(r[upwind_link_at_link == -1], 1.0)
    assert_array_almost_equal(
        r[upwind_link_at_link != -1], expected[upwind_link_at_link != -1]
    )

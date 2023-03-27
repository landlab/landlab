import numpy as np
from numpy.testing import assert_array_equal, assert_array_almost_equal

from landlab import HexModelGrid, RasterModelGrid
from landlab.components.advection import (
    find_upwind_link_at_link,
    upwind_to_local_grad_ratio,
)


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
    v = grid.add_zeros("value", at="node")
    v[:] = np.arange(grid.number_of_nodes) ** 2
    upwind_link_at_link = find_upwind_link_at_link(grid, 1.0)
    r = upwind_to_local_grad_ratio(grid, v, upwind_link_at_link)
    assert_array_almost_equal(
        r,
        [
            1.0,
            1.0,
            1.0,
            1.0,
            1.0,
            1.0,
            1.0,
            1.0,
            0.81818182,
            0.84615385,
            1.0,
            0.42857143,
            0.5,
            1.0,
            1.0,
            0.89473684,
            0.9047619,
            1.0,
            0.63636364,
            0.66666667,
            1.0,
            1.0,
            1.0,
            1.0,
        ],
    )

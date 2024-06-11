"""
Created on Sat Nov 14 10:48:23 2015

@author: gtucker
"""

import numpy as np
from numpy.testing import assert_array_equal

from landlab import RasterModelGrid


def test_link_order():
    """Test ordering of links."""
    rg = RasterModelGrid((3, 4))
    assert_array_equal(
        rg.node_at_link_tail, [0, 1, 2, 0, 1, 2, 3, 4, 5, 6, 4, 5, 6, 7, 8, 9, 10]
    )
    assert_array_equal(
        rg.node_at_link_head, [1, 2, 3, 4, 5, 6, 7, 5, 6, 7, 8, 9, 10, 11, 9, 10, 11]
    )


def test_link_at_face():
    grid = RasterModelGrid((4, 5), xy_spacing=1.0)
    assert_array_equal(
        grid.link_at_face,
        [5, 6, 7, 9, 10, 11, 12, 14, 15, 16, 18, 19, 20, 21, 23, 24, 25],
    )

    grid = RasterModelGrid((3, 4))
    assert_array_equal(grid.link_at_face, [4, 5, 7, 8, 9, 11, 12])


def test_horizontal_links():
    grid = RasterModelGrid((4, 5), xy_spacing=1.0)
    assert_array_equal(
        grid.horizontal_links,
        [0, 1, 2, 3, 9, 10, 11, 12, 18, 19, 20, 21, 27, 28, 29, 30],
    )


def test_vertical_links():
    grid = RasterModelGrid((4, 5), xy_spacing=1.0)
    assert_array_equal(
        grid.vertical_links, [4, 5, 6, 7, 8, 13, 14, 15, 16, 17, 22, 23, 24, 25, 26]
    )


def test_link_dirs_at_node():
    grid = RasterModelGrid((4, 3))
    assert_array_equal(
        grid.link_dirs_at_node,
        [
            [-1, -1, 0, 0],
            [-1, -1, 1, 0],
            [0, -1, 1, 0],
            [-1, -1, 0, 1],
            [-1, -1, 1, 1],
            [0, -1, 1, 1],
            [-1, -1, 0, 1],
            [-1, -1, 1, 1],
            [0, -1, 1, 1],
            [-1, 0, 0, 1],
            [-1, 0, 1, 1],
            [0, 0, 1, 1],
        ],
    )
    assert grid.link_dirs_at_node.dtype == np.int8

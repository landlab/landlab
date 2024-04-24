"""
Created on Sat Nov 14 10:36:03 2015

@author: gtucker
"""

import numpy as np
from numpy.testing import assert_array_equal
from pytest import approx

from landlab import HexModelGrid


def test_hex_grid_link_order():
    """Test the order of links in a small hex grid."""
    grid = HexModelGrid((3, 2))
    assert_array_equal(
        grid.nodes_at_link,
        [
            [0, 1],
            [0, 2],
            [0, 3],
            [1, 3],
            [1, 4],
            [2, 3],
            [3, 4],
            [2, 5],
            [3, 5],
            [3, 6],
            [4, 6],
            [5, 6],
        ],
    )

    grid = HexModelGrid((2, 3), orientation="vertical")
    assert_array_equal(
        grid.nodes_at_link,
        [
            [1, 0],
            [0, 2],
            [0, 3],
            [1, 3],
            [3, 2],
            [1, 4],
            [2, 5],
            [4, 3],
            [3, 5],
            [3, 6],
            [4, 6],
            [6, 5],
        ],
    )


def test_nodes_at_link():
    """Test nodes_at_link shares data with tail and head."""
    grid = HexModelGrid((3, 2))

    assert_array_equal(grid.nodes_at_link[:, 0], grid.node_at_link_tail)
    assert_array_equal(grid.nodes_at_link[:, 1], grid.node_at_link_head)

    assert np.may_share_memory(grid.nodes_at_link, grid.node_at_link_tail)
    assert np.may_share_memory(grid.nodes_at_link, grid.node_at_link_head)


def test_face_at_link():
    grid = HexModelGrid((3, 3))
    assert_array_equal(
        grid.face_at_link,
        [-1, -1, -1, 0, 1, 2, 3, -1, 4, 5, 6, -1, 7, 8, 9, 10, -1, -1, -1],
    )


def test_length_of_face():
    grid = HexModelGrid((3, 3))
    assert grid.length_of_face == approx(np.tan(np.pi / 6.0), abs=1e-5)
    assert len(grid.length_of_face) == grid.number_of_faces


def test_link_at_face():
    grid = HexModelGrid((3, 3))
    assert_array_equal(grid.link_at_face, [3, 4, 5, 6, 8, 9, 10, 12, 13, 14, 15])

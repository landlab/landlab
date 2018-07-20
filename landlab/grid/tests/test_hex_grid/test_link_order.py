# -*- coding: utf-8 -*-
"""
Created on Sat Nov 14 10:36:03 2015

@author: gtucker
"""
import numpy as np

from landlab import HexModelGrid
from numpy.testing import assert_array_equal


def test_hex_grid_link_order():
    """Test the order of links in a small hex grid."""
    hg = HexModelGrid(3, 2)
    assert_array_equal(hg.node_at_link_tail, [0, 0, 0, 1, 1, 2, 3, 2, 3, 3,
                                              4, 5])
    assert_array_equal(hg.node_at_link_head, [1, 2, 3, 3, 4, 3, 4, 5, 5, 6,
                                              6, 6])

    hg = HexModelGrid(2, 3, orientation='vertical')
    assert_array_equal(hg.node_at_link_tail, [1, 0, 0, 1, 3, 1, 2, 4, 3, 3,
                                              4, 6])
    assert_array_equal(hg.node_at_link_head, [0, 2, 3, 3, 2, 4, 5, 3, 5, 6,
                                              6, 5])


def test_nodes_at_link():
    """Test nodes_at_link shares data with tail and head."""
    grid = HexModelGrid(3, 2)

    assert_array_equal(grid.nodes_at_link[:, 0], grid.node_at_link_tail)
    assert_array_equal(grid.nodes_at_link[:, 1], grid.node_at_link_head)

    assert np.may_share_memory(grid.nodes_at_link, grid.node_at_link_tail)
    assert np.may_share_memory(grid.nodes_at_link, grid.node_at_link_head)


def test_face_at_link():
    grid = HexModelGrid(3, 3)
    assert_array_equal(grid.face_at_link,
                       [-1, -1, -1,  0,  1,  2,  3, -1,  4,  5,  6, -1, 7,  8,
                         9, 10, -1, -1, -1])


def test_width_of_face():
    grid = HexModelGrid(3, 3)
    assert_array_almost_equal(grid.width_of_face, np.tan(np.pi / 6.))
    assert_equal(len(grid.width_of_face), grid.number_of_faces)


def test_link_at_face():
    grid = HexModelGrid(3, 3)
    assert_array_equal(grid.link_at_face,
                       [ 3,  4,  5,  6,  8,  9, 10, 12, 13, 14, 15])

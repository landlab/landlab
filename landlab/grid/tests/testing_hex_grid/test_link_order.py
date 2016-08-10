# -*- coding: utf-8 -*-
"""
Created on Sat Nov 14 10:36:03 2015

@author: gtucker
"""

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


if __name__=='__main__':
    test_hex_grid_link_order()
    
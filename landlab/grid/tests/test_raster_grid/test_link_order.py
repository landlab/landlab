# -*- coding: utf-8 -*-
"""
Created on Sat Nov 14 10:48:23 2015

@author: gtucker
"""

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


if __name__ == "__main__":
    test_link_order()

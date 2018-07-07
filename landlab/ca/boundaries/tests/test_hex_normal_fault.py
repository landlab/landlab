#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Fri Jun 15 09:31:10 2018

@author: gtucker
"""

from landlab import HexModelGrid
from landlab.ca.boundaries.hex_lattice_tectonicizer import LatticeNormalFault
from numpy.testing import assert_array_equal


def test_links_to_update():
    """Test that update list includes lower 2 rows and fault-crossing links"""
    
    # Create a 6x6 test grid
    hg = HexModelGrid(6, 6, shape='rect', orientation='vert')

    lnf = LatticeNormalFault(grid=hg, fault_x_intercept=-0.1)
    
    assert_array_equal(lnf.links_to_update, [ 8,  9, 11, 12, 13, 14, 15, 16,
                                             18, 19, 20, 21, 22, 24, 25, 26,
                                             27, 28, 29, 30, 31, 32, 35, 36,
                                             37, 38, 40, 41, 43, 46, 51, 54,
                                             56, 60, 62, 68, 70, 73, 76, 77,
                                             80])

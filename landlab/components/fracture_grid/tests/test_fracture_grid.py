# -*- coding: utf-8 -*-
"""
test_fracture_grid.py

Unit test for fracture_grid component.

Created on Sun Oct 18 09:47:59 2015

@author: gtucker
"""

from landlab.components.fracture_grid import make_frac_grid
from numpy.testing import assert_array_equal


def test_frac_grid():
    
    frac_spacing = 3
    nrows = 9
    ncols = 9
    
    fg = make_frac_grid(frac_spacing, nrows, ncols, seed=1)
    
    assert_array_equal(fg, [[1, 0, 0, 0, 0, 0, 1, 1, 0],
                            [0, 0, 0, 0, 1, 1, 0, 0, 0],
                            [1, 1, 1, 1, 1, 1, 1, 0, 0],
                            [1, 1, 1, 0, 0, 0, 0, 1, 1],
                            [0, 0, 1, 1, 1, 1, 1, 1, 1],
                            [1, 1, 1, 1, 1, 0, 0, 0, 0],
                            [0, 0, 0, 0, 1, 1, 1, 0, 0],
                            [0, 0, 0, 0, 0, 1, 1, 1, 0],
                            [0, 0, 0, 0, 0, 0, 0, 1, 1]])
    
if __name__=='__main__':
    test_frac_grid()

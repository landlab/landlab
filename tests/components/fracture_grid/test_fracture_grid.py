# -*- coding: utf-8 -*-
"""
test_fracture_grid.py

Unit test for fracture_grid component.

Created on Sun Oct 18 09:47:59 2015

@author: gtucker
"""

from numpy.testing import assert_array_equal

from landlab import RasterModelGrid
from landlab.components.fracture_grid import FractureGridGenerator


def test_frac_grid():
    nrows = 9
    ncols = 9

    grid = RasterModelGrid((nrows, ncols))
    fg = FractureGridGenerator(grid, frac_spacing=3, seed=1)
    fg.run_one_step()

    assert_array_equal(
        grid.at_node["fracture_at_node"].reshape((9, 9)),
        [
            [1, 0, 0, 0, 0, 0, 1, 1, 0],
            [0, 0, 0, 0, 1, 1, 0, 0, 0],
            [1, 1, 1, 1, 1, 1, 1, 0, 0],
            [1, 1, 1, 0, 0, 0, 0, 1, 1],
            [0, 0, 1, 1, 1, 1, 1, 1, 1],
            [1, 1, 1, 1, 1, 0, 0, 0, 0],
            [0, 0, 0, 0, 1, 1, 1, 0, 0],
            [0, 0, 0, 0, 0, 1, 1, 1, 0],
            [0, 0, 0, 0, 0, 0, 0, 1, 1],
        ],
    )

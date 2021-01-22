# -*- coding: utf-8 -*-
"""
Created on Wed Jan 13 15:34:08 2021

@author: laure
"""

import numpy as np
from landlab import HexModelGrid
from landlab import RasterModelGrid

from landlab.utils.window_statistic import calculate_window_statistic

def test_hex_grid():
    grid = HexModelGrid((5, 4),spacing=10)
    grid.status_at_node[grid.status_at_node==1] = grid.BC_NODE_IS_CLOSED
    z = grid.add_zeros("topographic__elevation", at="node")
    z += np.arange(len(z))
    relief = calculate_window_statistic(grid,'topographic__elevation',
                                        np.ptp,search_radius=15,
                                        calc_on_closed_nodes=False)

    check = [     np.nan,  np.nan,  np.nan,  np.nan,
               np.nan,     6.,     7.,     7.,    np.nan,
            np.nan,    11.,    12.,    12.,    11.,  np.nan,
               np.nan,     7.,     7.,     6.,    np.nan,
                  np.nan,  np.nan,  np.nan,  np.nan]
    
    np.testing.assert_equal(relief, check)

def test_calc_with_interior_closed_nodes():
    grid = RasterModelGrid((5, 5), xy_spacing=10.0)
    grid.status_at_node[[7,8,12,13]] = grid.BC_NODE_IS_CLOSED
    z = grid.add_zeros("topographic__elevation", at="node")
    z += np.arange(len(z))
    relief = calculate_window_statistic(grid,'topographic__elevation',
                                        np.ptp,search_radius=15,
                                        calc_on_closed_nodes=False)

    check = [  6.,     6.,     5.,     7.,     6.,
              11.,    11.,  np.nan,  np.nan,   11.,
              11.,    12.,  np.nan,  np.nan,   10.,
              11.,    12.,    12.,     10.,    10.,
               6.,     7.,     7.,     7.,     6.]
    
    np.testing.assert_equal(relief, check)

# -*- coding: utf-8 -*-
"""
Created on Wed Jan 13 15:34:08 2021

@author: laure
"""

import numpy as np
from landlab import HexModelGrid

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

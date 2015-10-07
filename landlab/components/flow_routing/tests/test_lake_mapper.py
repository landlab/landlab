# -*- coding: utf-8 -*-
"""
test_lake_mapper: 

Created on Sun Sep 27 09:52:50, 2015

@author: gtucker
"""

import landlab
from landlab import RasterModelGrid
from landlab.components.flow_routing.lake_mapper import DepressionFinderAndRouter
from numpy import sin, pi
import numpy as np  # for use of np.round
from numpy.testing import assert_array_equal

NUM_GRID_ROWS = 8
NUM_GRID_COLS = 8
PERIOD_X = 8.
PERIOD_Y = 4.


def create_test_grid():
    """
    Create a test grid and elevation field with sinusoidal depressions and
    hills.
    """
    # Create grid
    rmg = RasterModelGrid(NUM_GRID_ROWS, NUM_GRID_COLS)
    
    # Create topography field
    z = rmg.add_zeros('node', 'topographic__elevation')
    
    # Make topography into sinusoidal hills and depressions
    z[:] = sin(2*pi*rmg.node_x/PERIOD_X) * sin(2*pi*rmg.node_y/PERIOD_Y)
    
    # Set 3 sides of the grid to be closed boundaries
    rmg.set_closed_boundaries_at_grid_edges(False, True, True, True)
    
    return rmg
    
    
def check_fields(grid):
    """
    Check to make sure the right fields have been created.
    """
    try:
        grid.at_node['topographic__elevation']
        grid.at_node['flood_status_code']
        grid.at_node['depression__depth']
        grid.at_node['depression__outlet_node_id']
        grid.at_node['is_pit']
    except:
        print('Test failure in check_fields')
        raise


def check_array_values(rmg, lm):
    """
    Check values of the various fields against known values.
    """
    assert_array_equal(lm.is_pit, \
    [False, False, False, False, False, False, False, False, False, False, False, False,
     False, False,  True, False, False, False, False, False, False, False, False, False,
     False, False,  True, False, False, False, False, False, False, False, False, False,
     False, False, False, False, False, False, False, False, False, False,  True, False,
     False,  True, False, False, False, False, False, False, False, False, False, False,
     False, False, False, False])

    assert_array_equal(lm.flood_status, \
    [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 3, 3, 0, 0, 0, 0, 0, 0, 3, 3, 0,
     0, 3, 3, 3, 3, 0, 0, 0, 0, 3, 3, 3, 3, 0, 0, 0, 0, 0, 0, 0, 0, 3, 3, 0,
     0, 3, 0, 0, 0, 3, 3, 0, 0, 0, 0, 0, 0, 0, 0, 0])

    dd1 = np.round(lm.depression_depth*100)
    assert_array_equal(dd1, \
    [ 0.,   0.,   0.,   0.,   0.,   0.,   0.,   0.,   0.,   0.,   0.,   0.,
      0.,  71., 100.,   0.,   0.,   0.,   0.,   0.,   0.,   0.,   0.,   0.,
      0.,  71., 100.,  71.,   0.,   0.,   0.,   0.,   0.,   0.,   0.,   0.,
      0.,   0.,   0.,   0.,   0.,   0.,   0.,   0.,   0.,  71., 100.,   0.,
      0.,   0.,   0.,   0.,   0.,   0.,   0.,   0.,   0.,   0.,   0.,   0.,
      0.,   0.,   0.,   0.])

    dd1 = np.round(rmg.at_node['depression__depth']*100)
    assert_array_equal(dd1, \
    [ 0.,   0.,   0.,   0.,   0.,   0.,   0.,   0.,   0.,   0.,   0.,   0.,
      0.,  71., 100.,   0.,   0.,   0.,   0.,   0.,   0.,   0.,   0.,   0.,
      0.,  71., 100.,  71.,   0.,   0.,   0.,   0.,   0.,   0.,   0.,   0.,
      0.,   0.,   0.,   0.,   0.,   0.,   0.,   0.,   0.,  71., 100.,   0.,
      0.,   0.,   0.,   0.,   0.,   0.,   0.,   0.,   0.,   0.,   0.,   0.,
      0.,   0.,   0.,   0.])
      
    assert_array_equal(lm.depression_outlet, \
    [2147483647, 2147483647, 2147483647, 2147483647, 2147483647, 2147483647,
     2147483647, 2147483647, 2147483647, 2147483647, 2147483647, 2147483647,
     2147483647,          5,          5, 2147483647, 2147483647, 2147483647,
     2147483647, 2147483647, 2147483647,          5,          5, 2147483647,
     2147483647,          5,          5,          5,          5, 2147483647,
     2147483647, 2147483647, 2147483647,          5,          5,          5,
              5, 2147483647, 2147483647, 2147483647, 2147483647, 2147483647,
     2147483647, 2147483647, 2147483647,          5,          5, 2147483647,
     2147483647,         50, 2147483647, 2147483647, 2147483647,          5,
              5, 2147483647, 2147483647, 2147483647, 2147483647, 2147483647,
     2147483647, 2147483647, 2147483647, 2147483647])

    assert_array_equal(rmg.at_node['depression__outlet_node_id'], \
    [2147483647, 2147483647, 2147483647, 2147483647, 2147483647, 2147483647,
     2147483647, 2147483647, 2147483647, 2147483647, 2147483647, 2147483647,
     2147483647,          5,          5, 2147483647, 2147483647, 2147483647,
     2147483647, 2147483647, 2147483647,          5,          5, 2147483647,
     2147483647,          5,          5,          5,          5, 2147483647,
     2147483647, 2147483647, 2147483647,          5,          5,          5,
              5, 2147483647, 2147483647, 2147483647, 2147483647, 2147483647,
     2147483647, 2147483647, 2147483647,          5,          5, 2147483647,
     2147483647,         50, 2147483647, 2147483647, 2147483647,          5,
              5, 2147483647, 2147483647, 2147483647, 2147483647, 2147483647,
     2147483647, 2147483647, 2147483647, 2147483647])


def test_lake_mapper():
    """
    Create a test grid and run a series of tests.
    """
    # Make a test grid
    rmg = create_test_grid()

    # Instantiate a lake mapper
    # (Note that we don't need to send it an input file name, because our grid
    # already has a topographic__elevation field)
    lm = DepressionFinderAndRouter(rmg)
    
    # Run it on our test grid
    lm.map_depressions()
    
    # Run tests
    check_fields(rmg)
    check_array_values(rmg, lm)
    
    
if __name__=='__main__':
    test_lake_mapper()
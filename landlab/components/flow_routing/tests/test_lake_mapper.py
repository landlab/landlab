# -*- coding: utf-8 -*-
"""
test_lake_mapper: 

Created on Sun Sep 27 09:52:50 2015

@author: gtucker
"""

from landlab import RasterModelGrid
from landlab.components.flow_routing.lake_mapper import DepressionFinderAndRouter
from numpy import sin, pi

NUM_GRID_ROWS = 32
NUM_GRID_COLS = 32
PERIOD_X = 32.
PERIOD_Y = 16.


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
    
#    for i in range(rmg.number_of_nodes/4):
#        print i, rmg.node_x[i], rmg.node_y[i], rmg.at_node['topographic__elevation'][i], \
#            lm.is_pit[i], lm.flood_status[i], rmg.at_node['depression__depth'][i], \
#            rmg.at_node['depression__outlet_node_id'][i]
            
    
    
if __name__=='__main__':
    test_lake_mapper()
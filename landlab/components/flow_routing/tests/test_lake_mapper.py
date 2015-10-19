# -*- coding: utf-8 -*-
"""
test_lake_mapper: 

Created on Sun Sep 27 09:52:50, 2015

@author: gtucker
"""

import landlab
from landlab import RasterModelGrid
from landlab.components.flow_routing.route_flow_dn import FlowRouter
from landlab.components.flow_routing.lake_mapper import DepressionFinderAndRouter
from numpy import sin, pi
import numpy as np  # for use of np.round
from numpy.testing import assert_array_equal
from landlab import BAD_INDEX_VALUE as XX

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

def setup_dans_grid():
    """
    Create a 7x7 test grid with a well defined hole in it.
    """
    from landlab import RasterModelGrid
    from landlab.components.flow_routing.route_flow_dn import FlowRouter
    from landlab.components.flow_routing.lake_mapper import DepressionFinderAndRouter
    
    mg = RasterModelGrid(7,7,1.)
    
    z = np.array([  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
                    0.0,  2.0,  2.0,  2.0,  2.0,  2.0,  0.0,
                    0.0,  2.0,  1.6,  1.5,  1.6,  2.0,  0.0,
                    0.0,  2.0,  1.7,  1.6,  1.7,  2.0,  0.0,
                    0.0,  2.0,  1.8,  2.0,  2.0,  2.0,  0.0,
                    0.0,  1.0,  0.6,  1.0,  1.0,  1.0,  0.0,
                    0.0,  0.0, -0.5,  0.0,  0.0,  0.0,  0.0])

    r_old = np.array([  0,  1,  2,  3,  4,  5,  6,
                        7,  1,  2,  3,  4,  5, 13,
                        14, 14, 17, 17, 17, 20, 20,
                        21, 21, 17, 17, 17, 27, 27,
                        28, 28, 37, 38, 39, 34, 34,
                        35, 44, 44, 44, 46, 47, 41,
                        42, 43, 44, 45, 46, 47, 48])

    r_new = np.array([  0,  1,  2,  3,  4,  5,  6,
                        7,  1,  2,  3,  4,  5, 13,
                        14, 14, 23, 23, 24, 20, 20,
                        21, 21, 30, 30, 24, 27, 27,
                        28, 28, 37, 38, 39, 34, 34,
                        35, 44, 44, 44, 46, 47, 41,
                        42, 43, 44, 45, 46, 47, 48])
    
    A_old = np.array([[  1.,  2.,  2.,  2.,  2.,  2.,  1.,
                         1.,  1.,  1.,  1.,  1.,  1.,  1.,
                         2.,  1.,  1.,  6.,  1.,  1.,  2.,
                         2.,  1.,  1.,  1.,  1.,  1.,  2.,
                         2.,  1.,  1.,  1.,  1.,  1.,  2.,
                         1.,  1.,  2.,  2.,  2.,  1.,  1.,
                         1.,  1.,  6.,  1.,  3.,  2.,  1.]])

    A_new = np.array([[  1.,  2.,  2.,  2.,  2.,  2.,  1.,
                         1.,  1.,  1.,  1.,  1.,  1.,  1.,
                         2.,  1.,  1.,  1.,  1.,  1.,  2.,
                         2.,  1.,  3.,  3.,  1.,  1.,  2.,
                         2.,  1.,  7.,  1.,  1.,  1.,  2.,
                         1.,  1.,  8.,  2.,  2.,  1.,  1.,
                         1.,  1., 12.,  1.,  3.,  2.,  1.]])
    
    s_new = np.array([  0,  1,  8,  2,  9,  3, 10,
                        4, 11,  5, 12,  6,  7, 13,
                       14, 15, 20, 19, 21, 22, 27,
                       26, 28, 29, 34, 33, 35, 41,
                       42, 43, 44, 36, 37, 30, 23,
                       16, 17, 24, 18, 25, 38, 31,
                       45, 46, 39, 32, 47, 40, 48])

    depression_outlet_target = np.array([ XX, XX, XX, XX, XX, XX, XX,
                                          XX, XX, XX, XX, XX, XX, XX,
                                          XX, XX, 30, 30, 30, XX, XX,
                                          XX, XX, 30, 30, 30, XX, XX,
                                          XX, XX, XX, XX, XX, XX, XX,
                                          XX, XX, XX, XX, XX, XX, XX,
                                          XX, XX, XX, XX, XX, XX, XX])
    
    mg.add_field('node', 'topographic__elevation', z, units='-')
    fr = FlowRouter(mg)
    lf = DepressionFinderAndRouter(mg)
    
    
def check_fields(grid):
    """
    Check to make sure the right fields have been created.
    """
    try:
        grid.at_node['topographic__elevation']
        grid.at_node['flood_status_code']
        grid.at_node['depression__depth']
        grid.at_node['depression__outlet_node']
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

    assert_array_equal(rmg.at_node['depression__outlet_node'], \
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

@with_setup(setup_dans_grid)
def test_initial_routing():
    """
    Test the action of fr.route_flow() on the grid.
    """
    fr.route_flow()
    assert_array_equal(mg.at_node['flow_receiver'], r_old)
    assert_array_equal(mg.at_node['drainage_area'], A_old)

@with_setup(setup_dans_grid)
def test_rerouting_with_supplied_pits():
    """
    Test with the output from a successful run of fr.route_flow.
    """
    fr.route_flow()
    lf.map_depressions()
    assert_array_equal(mg.at_node['flow_receiver'], r_new)
    assert_array_equal(mg.at_node['drainage_area'], A_new)
    assert_array_equal(mg.at_node['water__volume_flux'], A_new)  # as P=1
    assert_array_equal(mg.at_node['upstream_ID_order'], s_new)

@with_setup(setup_dans_grid)
def test_filling_alone():
    """
    Test the filler alone, w/o supplying information on the pits.
    """
    lf.map_depressions(pits=None, reroute_flow=False)
    assert_array_equal(mg.at_node['flow_receiver'], r_old)
    assert_array_equal(lf.depression_outlet, depression_outlet_target)

@with_setup(setup_dans_grid)
def test_filling_supplied_pits():
    """
    Test the filler without rereouting, but confusingly, where there *is*
    aready routing information available!
    Also tests the supply of an array for 'pits'
    """
    fr.route_flow()
    lf.map_depressions(pits=mg.at_node['flow_sinks'], reroute_flow=False)
    assert_array_equal(mg.at_node['flow_receiver'], r_old)

@with_setup(setup_dans_grid)
def test_pits_as_IDs():
    """
    Smoke test for passing specific IDs, not an array, to the mapper.
    """
    fr.route_flow()
    lf.map_depressions(pits=np.where(mg.at_node['flow_sinks'])[0])
    assert_array_equal(mg.at_node['drainage_area'], A_new)

def three_pits():
    """
    A test to ensure the component correctly handles cases where there are
    multiple pits.
    """
    mg = RasterModelGrid(10,10,1.)
    z = mg.add_field('node', 'topographic__elevation', mg.node_x.copy())
    np.random.seed(seed=0)
    z += np.random.rand(100)/10000.
    # punch some holes
    z[33] = 1.
    z[43] = 1.
    z[37] = 4.
    z[74:76] = 1.
    fr = FlowRouter(mg)
    lf = DepressionFinderAndRouter(mg)
    fr.route_flow()
    lf.map_depressions()
    
    
    
if __name__=='__main__':
    test_lake_mapper()
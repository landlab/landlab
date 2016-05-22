# -*- coding: utf-8 -*-
"""
test_lake_mapper:

Created on Sun Sep 27 09:52:50, 2015

@author: gtucker, amended dejh
"""

import landlab
from landlab import RasterModelGrid
from landlab.components.flow_routing import (FlowRouter,
                                             DepressionFinderAndRouter)
from numpy import sin, pi
import numpy as np  # for use of np.round
from numpy.testing import assert_array_equal, assert_array_almost_equal
from landlab import BAD_INDEX_VALUE as XX
from nose.tools import (with_setup, assert_true, assert_false,
                        assert_almost_equal, assert_equal)

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
    rmg.set_closed_boundaries_at_grid_edges(True, True, True, False)
    
    return rmg


def setup_dans_grid():
    """
    Create a 7x7 test grid with a well defined hole in it.
    """
    from landlab import RasterModelGrid
    from landlab.components.flow_routing import (FlowRouter,
                                                 DepressionFinderAndRouter)

    global fr, lf, mg
    global z, r_new, r_old, A_new, A_old, s_new, depr_outlet_target

    mg = RasterModelGrid(7, 7, 1.)

    z = np.array([0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
                  0.0,  2.0,  2.0,  2.0,  2.0,  2.0,  0.0,
                  0.0,  2.0,  1.6,  1.5,  1.6,  2.0,  0.0,
                  0.0,  2.0,  1.7,  1.6,  1.7,  2.0,  0.0,
                  0.0,  2.0,  1.8,  2.0,  2.0,  2.0,  0.0,
                  0.0,  1.0,  0.6,  1.0,  1.0,  1.0,  0.0,
                  0.0,  0.0, -0.5,  0.0,  0.0,  0.0,  0.0]).flatten()

    r_old = np.array([0,  1,  2,  3,  4,  5,  6,
                      7,  1,  2,  3,  4,  5, 13,
                     14, 14, 17, 17, 17, 20, 20,
                     21, 21, 17, 17, 17, 27, 27,
                     28, 28, 37, 38, 39, 34, 34,
                     35, 44, 44, 44, 46, 41, 41,
                     42, 43, 44, 45, 46, 47, 48]).flatten()

    r_new = np.array([0,  1,  2,  3,  4,  5,  6,
                      7,  1,  2,  3,  4,  5, 13,
                     14, 14, 23, 23, 24, 20, 20,
                     21, 21, 30, 30, 24, 27, 27,
                     28, 28, 37, 38, 39, 34, 34,
                     35, 44, 44, 44, 46, 41, 41,
                     42, 43, 44, 45, 46, 47, 48]).flatten()

    A_old = np.array([[0.,  1.,  1.,  1.,  1.,  1.,  0.,
                       0.,  1.,  1.,  1.,  1.,  1.,  0.,
                       1.,  1.,  1.,  6.,  1.,  1.,  1.,
                       1.,  1.,  1.,  1.,  1.,  1.,  1.,
                       1.,  1.,  1.,  1.,  1.,  1.,  1.,
                       0.,  1.,  2.,  2.,  2.,  1.,  1.,
                       0.,  0.,  5.,  0.,  2.,  0.,  0.]]).flatten()

    A_new = np.array([[0.,  1.,  1.,  1.,  1.,  1.,  0.,
                       0.,  1.,  1.,  1.,  1.,  1.,  0.,
                       1.,  1.,  1.,  1.,  1.,  1.,  1.,
                       1.,  1.,  3.,  3.,  1.,  1.,  1.,
                       1.,  1.,  7.,  1.,  1.,  1.,  1.,
                       0.,  1.,  8.,  2.,  2.,  1.,  1.,
                       0.,  0., 11.,  0.,  2.,  0.,  0.]]).flatten()

    s_new = np.array([0,  1,  8,  2,  9,  3, 10,
                      4, 11,  5, 12,  6,  7, 13,
                     14, 15, 20, 19, 21, 22, 27,
                     26, 28, 29, 34, 33, 35, 41,
                     40, 42, 43, 44, 36, 37, 30,
                     23, 16, 17, 24, 18, 25, 38,
                     31, 45, 46, 39, 32, 47, 48]).flatten()

    depr_outlet_target = np.array([XX, XX, XX, XX, XX, XX, XX,
                                   XX, XX, XX, XX, XX, XX, XX,
                                   XX, XX, 30, 30, 30, XX, XX,
                                   XX, XX, 30, 30, 30, XX, XX,
                                   XX, XX, XX, XX, XX, XX, XX,
                                   XX, XX, XX, XX, XX, XX, XX,
                                   XX, XX, XX, XX, XX, XX, XX]).flatten()

    mg.add_field('node', 'topographic__elevation', z, units='-')

    fr = FlowRouter(mg)
    lf = DepressionFinderAndRouter(mg)


def setup_D4_grid():
    """
    Test functionality of routing when D4 is specified.
    """
    global frD8, frD4, lfD8, lfD4, mg1, mg2
    global z, lake_nodes

    mg1 = RasterModelGrid(7, 7, 1.)
    mg2 = RasterModelGrid(7, 7, 1.)
    z = mg1.node_x.copy() + 1.
    lake_nodes = np.array([10, 16, 17, 18, 24, 32, 33, 38, 40])
    #z[lake_nodes] *= 0.01  #z[lake_nodes] = 0.
    z[lake_nodes] = 0.    
    mg1.add_field('node', 'topographic__elevation', z, units='-')
    mg2.add_field('node', 'topographic__elevation', z, units='-')

    frD8 = FlowRouter(mg1, method='D8')
    frD4 = FlowRouter(mg2, method='D4')
    lfD8 = DepressionFinderAndRouter(mg1, routing='D8')
    lfD4 = DepressionFinderAndRouter(mg2, routing='D4')

def check_fields1(grid):
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


def check_array_values1(rmg, lm):
    """
    Check values of the various fields against known values.
    """
#    for i in range(rmg.number_of_nodes):
#        print i, rmg.at_node['topographic__elevation'][i], lm.is_pit[i]

    assert_array_equal(lm.is_pit, \
    [False, False, False, False, False, False, False, False,
     False, False, False, False, False, False,  True, False,
     False, False, False, False, False, False, False, False,
     False, False,  True, False, False, False, False, False,
     False, False, False, False, False, False, False, False,
     False, False, False, False, False, False,  True, False,
     False,  True, False, False, False, False, False, False,
     False, False, False, False, False, False, False, False])

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
      
    assert_array_equal(lm.depression_outlet_map, \
    [XX, XX, XX, XX, XX, XX,
     XX, XX, XX, XX, XX, XX,
     XX,  5,  5, XX, XX, XX,
     XX, XX, XX,  5,  5, XX,
     XX,  5,  5,  5,  5, XX,
     XX, XX, XX,  5,  5,  5,
      5, XX, XX, XX, XX, XX,
     XX, XX, XX,  5,  5, XX,
     XX, 50, XX, XX, XX,  5,
      5, XX, XX, XX, XX, XX,
     XX, XX, XX, XX])

    assert_array_equal(rmg.at_node['depression__outlet_node'], \
    [XX, XX, XX, XX, XX, XX,
     XX, XX, XX, XX, XX, XX,
     XX,  5,  5, XX, XX, XX,
     XX, XX, XX,  5,  5, XX,
     XX,  5,  5,  5,  5, XX,
     XX, XX, XX,  5,  5,  5,
      5, XX, XX, XX, XX, XX,
     XX, XX, XX,  5,  5, XX,
     XX, 50, XX, XX, XX,  5,
      5, XX, XX, XX, XX, XX,
     XX, XX, XX, XX])


def setup_dans_grid2():
    """
    Create a 7x7 test grid with a well defined hole in it, AT THE EDGE.
    """
    from landlab import RasterModelGrid
    from landlab.components.flow_routing import (FlowRouter,
                                                 DepressionFinderAndRouter)

    global fr, lf, mg
    global z, r_new, r_old, A_new, A_old, s_new, depr_outlet_target

    mg = RasterModelGrid((7, 7), (1., 1.))

    z = mg.node_x.copy()
    guard_sides = np.concatenate((np.arange(7, 14), np.arange(35, 42)))
    edges = np.concatenate((np.arange(7), np.arange(42, 49)))
    hole_here = np.array(([15, 16, 22, 23, 29, 30]))
    z[guard_sides] = z[13]
    z[edges] = -2.  # force flow outwards from the tops of the guards
    z[hole_here] = -1.

    A_new = np.array([[[0.,   1.,   1.,   1.,   1.,   1.,   0.,
                        0.,   1.,   1.,   1.,   1.,   1.,   0.,
                       15.,   9.,   4.,   3.,   2.,   1.,   0.,
                        0.,   6.,   4.,   3.,   2.,   1.,   0.,
                        0.,   1.,   4.,   3.,   2.,   1.,   0.,
                        0.,   1.,   1.,   1.,   1.,   1.,   0.,
                        0.,   1.,   1.,   1.,   1.,   1.,   0.]]]).flatten()

    depr_outlet_target = np.array([XX, XX, XX, XX, XX, XX, XX,
                                   XX, XX, XX, XX, XX, XX, XX,
                                   XX, 14, 14, XX, XX, XX, XX,
                                   XX, 14, 14, XX, XX, XX, XX,
                                   XX, 14, 14, XX, XX, XX, XX,
                                   XX, XX, XX, XX, XX, XX, XX,
                                   XX, XX, XX, XX, XX, XX, XX]).flatten()

    mg.add_field('node', 'topographic__elevation', z, units='-')

    fr = FlowRouter(mg)
    lf = DepressionFinderAndRouter(mg)
    
def check_fields2(grid):
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


def check_array_values2(rmg, lm):
    """
    Check values of the various fields against known values.
    """
    assert_array_equal(lm.is_pit,
    [False, False, False, False, False, False, False, False, False, False, False, False,
     False, False,  True, False, False, False, False, False, False, False, False, False,
     False, False,  True, False, False, False, False, False, False, False, False, False,
     False, False, False, False, False, False, False, False, False, False,  True, False,
     False,  True, False, False, False, False, False, False, False, False, False, False,
     False, False, False, False])

    assert_array_equal(lm.flood_status,
    [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 3, 3, 0, 0, 0, 0, 0, 0, 3, 3, 0,
     0, 3, 3, 3, 3, 0, 0, 0, 0, 3, 3, 3, 3, 0, 0, 0, 0, 0, 0, 0, 0, 3, 3, 0,
     0, 3, 0, 0, 0, 3, 3, 0, 0, 0, 0, 0, 0, 0, 0, 0])

    dd1 = np.round(lm.depression_depth*100)
    assert_array_equal(dd1,
    [ 0.,   0.,   0.,   0.,   0.,   0.,   0.,   0.,   0.,   0.,   0.,   0.,
      0.,  71., 100.,   0.,   0.,   0.,   0.,   0.,   0.,   0.,   0.,   0.,
      0.,  71., 100.,  71.,   0.,   0.,   0.,   0.,   0.,   0.,   0.,   0.,
      0.,   0.,   0.,   0.,   0.,   0.,   0.,   0.,   0.,  71., 100.,   0.,
      0.,   0.,   0.,   0.,   0.,   0.,   0.,   0.,   0.,   0.,   0.,   0.,
      0.,   0.,   0.,   0.])

    dd1 = np.round(rmg.at_node['depression__depth']*100)
    assert_array_equal(dd1,
    [ 0.,   0.,   0.,   0.,   0.,   0.,   0.,   0.,   0.,   0.,   0.,   0.,
      0.,  71., 100.,   0.,   0.,   0.,   0.,   0.,   0.,   0.,   0.,   0.,
      0.,  71., 100.,  71.,   0.,   0.,   0.,   0.,   0.,   0.,   0.,   0.,
      0.,   0.,   0.,   0.,   0.,   0.,   0.,   0.,   0.,  71., 100.,   0.,
      0.,   0.,   0.,   0.,   0.,   0.,   0.,   0.,   0.,   0.,   0.,   0.,
      0.,   0.,   0.,   0.])
      
    assert_array_equal(lm.depression_outlet_map,
    [XX, XX, XX, XX, XX, XX,
     XX, XX, XX, XX, XX, XX,
     XX,  5,  5, XX, XX, XX,
     XX, XX, XX,  5,  5, XX,
     XX,  5,  5,  5,  5, XX,
     XX, XX, XX,  5,  5,  5,
      5, XX, XX, XX, XX, XX,
     XX, XX, XX,  5,  5, XX,
     XX, 50, XX, XX, XX,  5,
      5, XX, XX, XX, XX, XX,
     XX, XX, XX, XX])

    assert_array_equal(rmg.at_node['depression__outlet_node'], \
    [XX, XX, XX, XX, XX, XX,
     XX, XX, XX, XX, XX, XX,
     XX,  5,  5, XX, XX, XX,
     XX, XX, XX,  5,  5, XX,
     XX,  5,  5,  5,  5, XX,
     XX, XX, XX,  5,  5,  5,
      5, XX, XX, XX, XX, XX,
     XX, XX, XX,  5,  5, XX,
     XX, 50, XX, XX, XX,  5,
      5, XX, XX, XX, XX, XX,
     XX, XX, XX, XX])


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
    check_fields1(rmg)
    check_array_values1(rmg, lm)
    check_fields2(rmg)
    check_array_values2(rmg, lm)

@with_setup(setup_dans_grid)
def test_initial_routing():
    """
    Test the action of fr.route_flow() on the grid.
    """
    fr.route_flow()
    assert_array_equal(mg.at_node['flow__receiver_node'], r_old)
    assert_array_almost_equal(mg.at_node['drainage_area'], A_old)

@with_setup(setup_dans_grid)
def test_rerouting_with_supplied_pits():
    """
    Test with the output from a successful run of fr.route_flow.
    """
    fr.route_flow()
    lf.map_depressions()
    assert_array_equal(mg.at_node['flow__receiver_node'], r_new)
    assert_array_almost_equal(mg.at_node['drainage_area'], A_new)
    assert_array_almost_equal(mg.at_node['water__discharge'], A_new)
    assert_array_equal(mg.at_node['flow__upstream_node_order'], s_new)

@with_setup(setup_dans_grid)
def test_filling_alone():
    """
    Test the filler alone, w/o supplying information on the pits.
    """
    lf.map_depressions(pits=None, reroute_flow=False)
    assert_array_equal(mg.at_node['flow__receiver_node'], np.zeros(49, dtype=float))
    assert_array_equal(lf.depression_outlet_map, depr_outlet_target)

@with_setup(setup_dans_grid)
def test_filling_supplied_pits():
    """
    Test the filler without rereouting, but confusingly, where there *is*
    aready routing information available!
    Also tests the supply of an array for 'pits'
    """
    fr.route_flow()
    lf.map_depressions(pits=mg.at_node['flow__sink_flag'], reroute_flow=False)
    assert_array_equal(mg.at_node['flow__receiver_node'], r_old)

@with_setup(setup_dans_grid)
def test_pits_as_IDs():
    """
    Smoke test for passing specific IDs, not an array, to the mapper.
    """
    fr.route_flow()
    lf.map_depressions(pits=np.where(mg.at_node['flow__sink_flag'])[0])
    assert_array_almost_equal(mg.at_node['drainage_area'], A_new)


@with_setup(setup_dans_grid2)
def test_edge_draining():
    """
    This tests when the lake attempts to drain from an edge, where an issue
    is suspected.
    """
    fr.route_flow()
    lf.map_depressions()
    assert_array_almost_equal(mg.at_node['drainage_area'], A_new)
    assert_array_equal(lf.depression_outlet_map, depr_outlet_target)

def test_degenerate_drainage():
    """
    This "hourglass" configuration should be one of the hardest to correctly
    re-route.
    """
    mg = RasterModelGrid(9,5)
    z_init = mg.node_x.copy()*0.0001 + 1.
    lake_pits = np.array([7, 11, 12, 13, 17, 27, 31, 32, 33, 37])
    z_init[lake_pits] = -1.
    z_init[22] = 0.  # the common spill pt for both lakes
    z_init[21] = 0.1  # an adverse bump in the spillway
    z_init[20] = -0.2  # the spillway
    z = mg.add_field('node', 'topographic__elevation', z_init)

    fr = FlowRouter(mg)
    lf = DepressionFinderAndRouter(mg)
    fr.route_flow()
    lf.map_depressions()

    correct_A = np.array([ 0.,   0.,   0.,   0.,   0.,
                           0.,   1.,   3.,   1.,   0.,
                           0.,   5.,   1.,   2.,   0.,
                           0.,   1.,  10.,   1.,   0.,
                          21.,  21.,   1.,   1.,   0.,
                           0.,   1.,   9.,   1.,   0.,
                           0.,   3.,   1.,   2.,   0.,
                           0.,   1.,   1.,   1.,   0.,
                           0.,   0.,   0.,   0.,   0.])
    
    thelake = np.concatenate((lake_pits, [22])).sort()

    assert_array_almost_equal(mg.at_node['drainage_area'], correct_A)
    
    # assert np.all(np.equal(lf.lake_map[thelake], lf.lake_map[thelake[0]]))
    # assert not lf.lake_map[thelake[0]] == XX

def test_three_pits():
    """
    A test to ensure the component correctly handles cases where there are
    multiple pits.
    """
    mg = RasterModelGrid(10,10,1.)
    z = mg.add_field('node', 'topographic__elevation', mg.node_x.copy())
    # a sloping plane
    #np.random.seed(seed=0)
    #z += np.random.rand(100)/10000.
    # punch some holes
    z[33] = 1.
    z[43] = 1.
    z[37] = 4.
    z[74:76] = 1.
    fr = FlowRouter(mg)
    lf = DepressionFinderAndRouter(mg)
    fr.route_flow()
    lf.map_depressions()
    
    flow_sinks_target = np.zeros(100, dtype=bool)
    flow_sinks_target[mg.boundary_nodes] = True
    # no internal sinks now:
    assert_array_equal(mg.at_node['flow__sink_flag'], flow_sinks_target)
    
    # test conservation of mass:
    assert_almost_equal(mg.at_node['drainage_area'
                                       ].reshape((10,10))[1:-1,1].sum(), 8.**2)
    # ^all the core nodes
    
    # test the actual flow field:
    nA = np.array([  0.,   0.,   0.,   0.,   0.,   0.,   0.,   0.,   0.,   0.,
                     8.,   8.,   7.,   6.,   5.,   4.,   3.,   2.,   1.,   0.,
                     2.,   2.,   1.,   1.,   2.,   1.,   1.,   1.,   1.,   0.,
                    26.,  26.,  25.,  15.,  11.,  10.,   9.,   8.,   1.,   0.,
                     2.,   2.,   1.,   9.,   2.,   1.,   1.,   1.,   1.,   0.,
                     2.,   2.,   1.,   1.,   5.,   4.,   3.,   2.,   1.,   0.,
                     2.,   2.,   1.,   1.,   1.,   1.,   3.,   2.,   1.,   0.,
                    20.,  20.,  19.,  18.,  17.,  12.,   3.,   2.,   1.,   0.,
                     2.,   2.,   1.,   1.,   1.,   1.,   3.,   2.,   1.,   0.,
                     0.,   0.,   0.,   0.,   0.,   0.,   0.,   0.,   0.,   0.])
    assert_array_equal(mg.at_node['drainage_area'], nA)
    
    #test a couple more properties:
    lc = np.empty(100, dtype=int)
    lc.fill(XX)
    lc[33] = 33
    lc[43] = 33
    lc[37] = 37
    lc[74:76] = 74
    assert_array_equal(lf.lake_map, lc)
    assert_array_equal(lf.lake_codes, [33, 37, 74])
    assert_equal(lf.number_of_lakes, 3)
    assert_array_almost_equal(lf.lake_areas, [2., 1., 2.])
    assert_array_almost_equal(lf.lake_volumes, [2., 2., 4.])

def test_composite_pits():
    """
    A test to ensure the component correctly handles cases where there are
    multiple pits, inset into each other.
    """
    mg = RasterModelGrid(10, 10, 1.)
    z = mg.add_field('node', 'topographic__elevation', mg.node_x.copy())
    # a sloping plane
    #np.random.seed(seed=0)
    #z += np.random.rand(100)/10000.
    # punch one big hole
    z.reshape((10,10))[3:8,3:8] = 0.
    # dig a couple of inset holes
    z[57] = -1.
    z[44] = -2.
    z[54] = -10.
    fr = FlowRouter(mg)
    lf = DepressionFinderAndRouter(mg)
    fr.route_flow()
    lf.map_depressions()
    
    flow_sinks_target = np.zeros(100, dtype=bool)
    flow_sinks_target[mg.boundary_nodes] = True
    # no internal sinks now:
    assert_array_equal(mg.at_node['flow__sink_flag'], flow_sinks_target)
    
    # test conservation of mass:
    assert_almost_equal(mg.at_node['drainage_area'
                                       ].reshape((10,10))[1:-1,1].sum(), 8.**2)
    # ^all the core nodes
    
    # test the actual flow field:
    nA = np.array([  0.,   0.,   0.,   0.,   0.,   0.,   0.,   0.,   0.,   0.,
                     8.,   8.,   7.,   6.,   5.,   4.,   3.,   2.,   1.,   0.,
                     1.,   1.,   1.,   1.,   1.,   1.,   1.,   1.,   1.,   0.,
                     1.,   1.,   1.,   4.,   2.,   2.,   8.,   4.,   1.,   0.,
                     1.,   1.,   1.,   8.,   3.,  15.,   3.,   2.,   1.,   0.,
                     1.,   1.,   1.,  13.,  25.,   6.,   3.,   2.,   1.,   0.,
                     1.,   1.,   1.,  45.,   3.,   3.,   5.,   2.,   1.,   0.,
                    50.,  50.,  49.,   3.,   2.,   2.,   2.,   4.,   1.,   0.,
                     1.,   1.,   1.,   1.,   1.,   1.,   1.,   1.,   1.,   0.,
                     0.,   0.,   0.,   0.,   0.,   0.,   0.,   0.,   0.,   0.])
    assert_array_equal(mg.at_node['drainage_area'], nA)
    
    # the lake code map:
    lc = np.array([XX, XX, XX, XX, XX, XX, XX, XX, XX, XX,
                   XX, XX, XX, XX, XX, XX, XX, XX, XX, XX,
                   XX, XX, XX, XX, XX, XX, XX, XX, XX, XX,
                   XX, XX, XX, 57, 57, 57, 57, 57, XX, XX,
                   XX, XX, XX, 57, 57, 57, 57, 57, XX, XX,
                   XX, XX, XX, 57, 57, 57, 57, 57, XX, XX,
                   XX, XX, XX, 57, 57, 57, 57, 57, XX, XX,
                   XX, XX, XX, 57, 57, 57, 57, 57, XX, XX,
                   XX, XX, XX, XX, XX, XX, XX, XX, XX, XX,
                   XX, XX, XX, XX, XX, XX, XX, XX, XX, XX])
    
    #test the remaining properties:
    assert_equal(lf.lake_outlets.size, 1)
    assert_equal(lf.lake_outlets[0], 72)
    outlets_in_map = np.unique(lf.depression_outlet_map)
    assert_equal(outlets_in_map.size, 2)
    assert_equal(outlets_in_map[1], 72)
    assert_equal(lf.number_of_lakes, 1)
    assert_equal(lf.lake_codes[0], 57)
    assert_array_equal(lf.lake_map, lc)
    assert_almost_equal(lf.lake_areas[0], 25.)
    assert_almost_equal(lf.lake_volumes[0], 63.)

@with_setup(setup_D4_grid)
def test_D8_D4_fill():
    """
    Tests the functionality of D4 filling.
    """
    lfD8.map_depressions(pits=None, reroute_flow=False)
    lfD4.map_depressions(pits=None, reroute_flow=False)
    assert_equal(lfD8.number_of_lakes, 1)
    assert_equal(lfD4.number_of_lakes, 3)
    
    correct_D8_lake_map = np.empty(7*7, dtype=int)
    correct_D8_lake_map.fill(XX)
    correct_D8_lake_map[lake_nodes] = 10
    correct_D4_lake_map = correct_D8_lake_map.copy()
    correct_D4_lake_map[lake_nodes[5:]] = 32
    correct_D4_lake_map[lake_nodes[-2]] = 38
    correct_D8_depths = np.zeros(7*7, dtype=float)
    correct_D8_depths[lake_nodes] = 2.
    correct_D4_depths = correct_D8_depths.copy()
    correct_D4_depths[lake_nodes[5:]] = 4.
    correct_D4_depths[lake_nodes[-2]] = 3.
    
    assert_array_equal(lfD8.lake_map, correct_D8_lake_map)
    assert_array_equal(lfD4.lake_map, correct_D4_lake_map)
    
    assert_array_almost_equal(mg1.at_node['depression__depth'],
                              correct_D8_depths)
    assert_array_almost_equal(mg2.at_node['depression__depth'],
                              correct_D4_depths)

@with_setup(setup_D4_grid)
def test_D8_D4_route():
    """
    Tests the functionality of D4 routing.
    """
    frD8.route_flow()
    frD4.route_flow()
    lfD8.map_depressions()
    lfD4.map_depressions()
    assert_equal(lfD8.number_of_lakes, 1)
    assert_equal(lfD4.number_of_lakes, 3)

    flow_recD8 = np.array([ 0,  1,  2,  3,  4,  5,  6,  7, 16, 10, 16, 10, 18,
                           13, 14, 14, 15, 16, 10, 18, 20, 21, 16, 16, 16, 18,
                           33, 27, 28, 28, 24, 24, 24, 32, 34, 35, 35, 38, 32,
                           32, 32, 41, 42, 43, 44, 45, 46, 47, 48])
    flow_recD4 = np.array([ 0,  1,  2,  3,  4,  5,  6,  7,  7, 10, 17, 10, 11,
                           13, 14, 14, 15, 16, 17, 18, 20, 21, 21, 16, 17, 18,
                           33, 27, 28, 28, 29, 24, 31, 32, 34, 35, 35, 36, 37,
                           32, 33, 41, 42, 43, 44, 45, 46, 47, 48])
    assert_array_equal(mg1.at_node['flow__receiver_node'], flow_recD8)
    assert_array_equal(mg2.at_node['flow__receiver_node'], flow_recD4)
    assert_array_almost_equal(mg1.at_node['drainage_area'].reshape((7,7))[:,
                                  0].sum(),
                              mg2.at_node['drainage_area'].reshape((7,7))[:,
                                  0].sum())


if __name__=='__main__':
    setup_dans_grid()
    test_initial_routing()
    test_pits_as_IDs()
    test_rerouting_with_supplied_pits()

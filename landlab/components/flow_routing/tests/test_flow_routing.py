"""Test the flow router component.

This just tests the functionality of the router in toto - no attempt is made
to test the submodules.
Sinks are tested as part of the lake_mapper.

@author: dejh
"""
# Created on Thurs Nov 12, 2015
import os

from six.moves import range

import numpy as np
from numpy.testing import assert_array_equal, assert_array_almost_equal

from nose.tools import (with_setup, assert_true, assert_false, assert_raises,
                        assert_almost_equal, assert_equal)
try:
    from nose.tools import (assert_is, assert_set_equal, assert_dict_equal)
except ImportError:
    from landlab.testing.tools import (assert_is, assert_set_equal,
                                       assert_dict_equal)

import landlab
from landlab import RasterModelGrid, RadialModelGrid, FieldError
from landlab.components.flow_routing import FlowRouter
from landlab import CLOSED_BOUNDARY
from landlab import BAD_INDEX_VALUE as XX


_THIS_DIR = os.path.abspath(os.path.dirname(__file__))


def setup_dans_grid1():
    """
    Create a 5x5 test grid.
    This is a sheet flow test.
    """
    global fr, mg, infile
    global z, Q_in, A_target, frcvr_target, upids_target, Q_target, \
        steepest_target, links2rcvr_target

    mg = RasterModelGrid((5, 5), spacing=(10., 10.))

    infile = os.path.join(_THIS_DIR, 'test_fr_input.txt')

    z = mg.node_x.copy()

    Q_in = np.full(25, 2.)

    A_target = np.array([0.,  0.,  0.,  0.,  0.,
                         3.,  3.,  2.,  1.,  0.,
                         3.,  3.,  2.,  1.,  0.,
                         3.,  3.,  2.,  1.,  0.,
                         0.,  0.,  0.,  0.,  0.])*100.

    frcvr_target = np.array([0,  1,  2,  3,  4,
                             5,  5,  6,  7,  9,
                            10, 10, 11, 12, 14,
                            15, 15, 16, 17, 19,
                            20, 21, 22, 23, 24])

    upids_target = np.array([0,  1,  2,  3,  4,
                             5,  6,  7,  8,  9,
                            10, 11, 12, 13, 14,
                            15, 16, 17, 18, 19,
                            20, 21, 22, 23, 24])

    links2rcvr_target = np.full(25, XX)
    links2rcvr_target[mg.core_nodes] = np.array([ 9, 10, 11,
                                                 18, 19, 20,
                                                 27, 28, 29])

    Q_target = A_target * 2.  # only once Q_in is used

    steepest_target = np.array([0.,  0.,  0.,  0.,  0.,
                                0.,  1.,  1.,  1.,  0.,
                                0.,  1.,  1.,  1.,  0.,
                                0.,  1.,  1.,  1.,  0.,
                                0.,  0.,  0.,  0.,  0.])

    mg.add_field('node', 'topographic__elevation', z, units='-')


def setup_dans_grid2():
    """
    Create a 5x5 test grid.
    This tests more complex routing, with diffs between D4 & D8.
    """
    global fr, mg, infile
    global z, A_target_D8, A_target_D4, frcvr_target_D8, frcvr_target_D4, \
        upids_target_D8, upids_target_D4, steepest_target_D8, \
        steepest_target_D4, links2rcvr_target_D8, links2rcvr_target_D4

    mg = RasterModelGrid((5, 5), spacing=(10., 10.))

    infile = os.path.join(_THIS_DIR, 'test_fr_input.txt')

    z = np.array([7.,  7.,  7.,  7.,  7.,
                  7.,  5., 3.2,  6.,  7.,
                  7.,  2.,  3.,  5.,  7.,
                  7.,  1., 1.9,  4.,  7.,
                  7.,  0.,  7.,  7.,  7.])

    A_target_D8 = np.array([0.,     0.,     0.,     0.,     0.,
                            0.,   100.,   200.,   100.,     0.,
                            0.,   400.,   100.,   100.,     0.,
                            0.,   600.,   300.,   100.,     0.,
                            0.,   900.,     0.,     0.,     0.])

    A_target_D4 = np.array([0.,     0.,     0.,     0.,     0.,
                            0.,   100.,   200.,   100.,     0.,
                            0.,   200.,   400.,   100.,     0.,
                            0.,   900.,   600.,   100.,     0.,
                            0.,   900.,     0.,     0.,     0.])

    frcvr_target_D8 = np.array([0, 1, 2, 3, 4, 5, 11, 11, 7, 9, 10, 16, 16, 17,
                                14, 15, 21, 21, 17, 19, 20, 21, 22, 23, 24])

    frcvr_target_D4 = np.array([ 0,  1,  2,  3,  4, 
                                 5, 11, 12,  7,  9,
                                10, 16, 17, 12, 14,
                                15, 21, 16, 17, 19,
                                20, 21, 22, 23, 24])

    upids_target_D8 = np.array([0, 1, 2, 3, 4, 5, 9, 10, 14, 15, 19, 20, 21,
                                16, 11, 6, 7, 8, 12, 17, 13, 18, 22, 23, 24])

    upids_target_D4 = np.array([0, 1, 2, 3, 4, 5, 9, 10, 14, 15, 19, 20, 21,
                                16, 11, 6, 17, 12, 7, 8, 13, 18, 22, 23, 24])

    links2rcvr_target_D8 = np.full(25, XX)
    links2rcvr_target_D8[mg.core_nodes] = np.array([14, 51, 11,
                                                    23, 59, 61,
                                                    32, 67, 29])

    links2rcvr_target_D4 = np.full(25, XX)
    links2rcvr_target_D4[mg.core_nodes] = np.array([14, 15, 11,
                                                    23, 24, 20,
                                                    32, 28, 29])

    steepest_target_D8 = np.array([0., 0., 0., 0., 0.,
                                   0., 0.3, 0.08485281, 0.28, 0.,
                                   0., 0.1, 0.14142136, 0.21920310, 0.,
                                   0., 0.1, 0.13435029, 0.21,  0.,
                                   0., 0., 0., 0.,  0.])

    steepest_target_D4 = np.array([0., 0., 0., 0., 0.,
                                   0., 0.3, 0.02, 0.28, 0.,
                                   0., 0.1, 0.11, 0.2, 0.,
                                   0., 0.1, 0.09, 0.21, 0.,
                                   0., 0., 0., 0., 0.])

    mg.add_field('node', 'topographic__elevation', z, units='-')


def setup_internal_closed():
    """
    Create a 6x5 test grid, but with two internal nodes closed.
    This is a sheet flow test.
    """
    global fr, mg
    global z, Q_in, A_target, frcvr_target, upids_target, Q_target, \
        steepest_target, links2rcvr_target

    mg = RasterModelGrid((6, 5), spacing=(10., 10.))

    mg.set_closed_boundaries_at_grid_edges(True, True, False, True)
    mg.status_at_node[7] = CLOSED_BOUNDARY
    mg.status_at_node[16] = CLOSED_BOUNDARY

    z = mg.node_x.copy()

    Q_in = np.full(25, 2.)

    A_target = np.array([0.,  0.,  0.,  0.,  0.,
                         1.,  1.,  0.,  1.,  0.,
                         6.,  6.,  3.,  1.,  0.,
                         0.,  0.,  2.,  1.,  0.,
                         3.,  3.,  2.,  1.,  0.,
                         0.,  0.,  0.,  0.,  0.])*100.

    frcvr_target = np.array([0,  1,  2,  3,  4,
                             5,  5,  7, 12,  9,
                            10, 10, 11, 12, 14,
                            15, 16, 11, 17, 19,
                            20, 20, 21, 22, 24,
                            25, 26, 27, 28, 29])

    links2rcvr_target = np.full(mg.number_of_nodes, XX)
    links2rcvr_target[mg.core_nodes] = np.array([ 9, 62,
                                                 18, 19, 20,
                                                 67, 29,
                                                 36, 37, 38])

    steepest_target = np.array([0.,  0.,  0.,  0.,  0.,
                                0.,  1.,  0.,  0.,  0.,
                                0.,  1.,  1.,  1.,  0.,
                                0.,  0.,  0.,  1.,  0.,
                                0.,  1.,  1.,  1.,  0.,
                                0.,  0.,  0.,  0.,  0.])
    steepest_target[np.array([8, 17])] = 1./np.sqrt(2.)

    mg.add_field('node', 'topographic__elevation', z, units='-')


def setup_voronoi():
    """
    Setup a simple 20 point Voronoi Delaunay grid (radial for ease)
    """
    global fr, vmg
    global A_target_core, A_target_outlet
    vmg = RadialModelGrid(2, dr=2.)
    z = np.full(20, 10., dtype=float)
    #vmg.status_at_node[8:] = CLOSED_BOUNDARY
    all_bounds_but_one = np.array((0, 1, 2, 3, 4, 7, 11, 15, 16, 17, 18, 19))
    vmg.status_at_node[all_bounds_but_one] = CLOSED_BOUNDARY
    #z[7] = 0.  # outlet
    z[12] = 0.  # outlet
    #inner_elevs = (3., 1., 4., 5., 6., 7., 8.)
    inner_elevs = (8., 7., 3., 1., 6., 4., 5.)
    #z[:7] = np.array(inner_elevs)
    z[vmg.core_nodes] = np.array(inner_elevs)
    vmg.add_field('node', 'topographic__elevation', z, units='-')
    fr = FlowRouter(vmg)

#    nodes_contributing = [np.array([0, 3, 4, 5]),
#                          np.array([0, 1, 2, 3, 4, 5, 6]),
#                          np.array([2, ]),
#                          np.array([3, ]),
#                          np.array([4, ]),
#                          np.array([5, ]),
#                          np.array([6, ])]

    # The follow list contains arrays with the IDs of cells contributing flow
    # to nodes 5, 6, 8, 9, 10, 13, and 14, respectively (which correspond to
    # cells 0-6)
    cells_contributing = [np.array([0, ]),
                          np.array([1, ]),
                          np.array([1, 2, 4, 6]),
                          np.array([0, 1, 2, 3, 4, 5, 6]),
                          np.array([4, ]),
                          np.array([5, ]),
                          np.array([6, ])]

    A_target_core = np.zeros(vmg.number_of_core_nodes)
    for i in range(7):
        A_target_core[i] = vmg.area_of_cell[cells_contributing[i]].sum()
    A_target_outlet = vmg.area_of_cell.sum()


def setup_voronoi_closedinternal():
    """
    Close one of the formerly core internal nodes.
    """
    global fr, vmg
    global A_target_internal, A_target_outlet
    vmg = RadialModelGrid(2, dr=2.)
    z = np.full(20, 10., dtype=float)
#    vmg.status_at_node[8:] = CLOSED_BOUNDARY
#    vmg.status_at_node[0] = CLOSED_BOUNDARY  # new internal closed
#    z[7] = 0.  # outlet
#    inner_elevs = (3., 1., 4., 5., 6., 7., 8.)
#    z[:7] = np.array(inner_elevs)
    all_bounds_but_one = np.array((0, 1, 2, 3, 4, 7, 11, 15, 16, 17, 18, 19))
    vmg.status_at_node[all_bounds_but_one] = CLOSED_BOUNDARY
    vmg.status_at_node[8] = CLOSED_BOUNDARY  # new internal closed
    z[12] = 0.  # outlet
    inner_elevs = (8., 7., 1., 6., 4., 5.)
    z[vmg.core_nodes] = np.array(inner_elevs)
    vmg.add_field('node', 'topographic__elevation', z, units='-')
    fr = FlowRouter(vmg)

#    nodes_contributing = [[],
#                          np.array([1, 2, 3, 4, 5, 6]),
#                          np.array([2, 3, 4, 5]),
#                          np.array([3, 4, 5]),
#                          np.array([4, 5]),
#                          np.array([5, ]),
#                          np.array([6, ])]

    cells_contributing = [np.array([0, ]),
                          np.array([1, ]),
                          np.array([0, 1, 3, 4, 5, 6]),
                          np.array([1, 4]),
                          np.array([1, 4, 5, 6]),
                          np.array([1, 4, 6])]

    A_target_internal = np.zeros(vmg.number_of_core_nodes, dtype=float)
    for i in range(6):
        A_target_internal[i] = vmg.area_of_cell[cells_contributing[i]].sum()
#        A_target_internal[i] = vmg.area_of_cell[nodes_contributing[i]].sum()
    A_target_outlet = vmg.area_of_cell[vmg.cell_at_node[vmg.core_nodes]].sum()


@with_setup(setup_dans_grid1)
def test_check_fields():
    """Check to make sure the right fields have been created."""
    fr = FlowRouter(mg)
    assert_array_equal(z, mg.at_node['topographic__elevation'])
    assert_array_equal(np.zeros(25), mg.at_node['drainage_area'])
    assert_array_equal(np.ones(25), mg.at_node['water__unit_flux_in'])

    fr = FlowRouter(mg, infile)
    assert_array_equal(np.full(25, 2.), mg.at_node['water__unit_flux_in'])


@with_setup(setup_dans_grid1)
def test_check_field_input():
    """Check we can successfully pass water__discharge_in."""
    mg.add_field('node', 'water__unit_flux_in',
                 np.full(25, 3.), units='m**3/s')
    fr = FlowRouter(mg)
    assert_array_equal(np.full(25, 3.), mg.at_node['water__unit_flux_in'])
    fr = FlowRouter(mg, infile)
    assert_array_equal(np.full(25, 2.), mg.at_node['water__unit_flux_in'])


@with_setup(setup_dans_grid1)
def test_accumulate_D8():
    """Test accumulation works for D8 in a simple scenario."""
    fr = FlowRouter(mg)
    fr.route_flow()
    assert_array_equal(A_target, mg.at_node['drainage_area'])
    assert_array_equal(frcvr_target, mg.at_node['flow__receiver_node'])
    assert_array_equal(upids_target, mg.at_node['flow__upstream_node_order'])
    assert_array_equal(links2rcvr_target, mg.at_node['flow__link_to_receiver_node'])
    assert_array_equal(A_target, mg.at_node['water__discharge'])
    assert_array_equal(steepest_target,
                       mg.at_node['topographic__steepest_slope'])


@with_setup(setup_dans_grid1)
def test_variable_Qin():
    """Test variable Qin field."""
    Qin_local = np.zeros(25, dtype=float)
    Qin_local[13] = 2.
    mg.add_field('node', 'water__unit_flux_in',
                 Qin_local, units='m**3/s')
    fr = FlowRouter(mg)
    fr.route_flow()
    Qout_local = np.zeros_like(Qin_local)
    Qout_local[10:14] = 200.
    assert_array_equal(Qout_local, mg.at_node['water__discharge'])
    assert_array_equal(A_target, mg.at_node['drainage_area'])
    # note that A DOES NOT CHANGE when messing with Q_in


@with_setup(setup_dans_grid2)
def test_irreg_topo():
    """Test D8 routing on a toy irregular topo."""
    fr = FlowRouter(mg)
    fr.route_flow()
    assert_array_equal(A_target_D8, mg.at_node['drainage_area'])
    assert_array_equal(frcvr_target_D8, mg.at_node['flow__receiver_node'])
    assert_array_equal(upids_target_D8, mg.at_node['flow__upstream_node_order'])
    assert_array_equal(links2rcvr_target_D8,
                       mg.at_node['flow__link_to_receiver_node'])
    assert_array_almost_equal(steepest_target_D8,
                              mg.at_node['topographic__steepest_slope'])


@with_setup(setup_dans_grid2)
def test_irreg_topo_old():
    """
    Test D4 routing on a toy irregular topo, old style, where 'method' is
    passed to the run method, not the init.
    """
    fr = FlowRouter(mg)
    fr.route_flow(method='D4')
    assert_array_equal(A_target_D4, mg.at_node['drainage_area'])
    assert_array_equal(frcvr_target_D4, mg.at_node['flow__receiver_node'])
    assert_array_equal(upids_target_D4, mg.at_node['flow__upstream_node_order'])
    assert_array_equal(links2rcvr_target_D4,
                       mg.at_node['flow__link_to_receiver_node'])
    assert_array_almost_equal(steepest_target_D4,
                              mg.at_node['topographic__steepest_slope'])


@with_setup(setup_dans_grid2)
def test_irreg_topo_new():
    """Test D4 routing on a toy irregular topo. 'method' passed to init."""
    fr = FlowRouter(mg, method='D4')
    fr.route_flow()
    assert_array_equal(A_target_D4, mg.at_node['drainage_area'])
    assert_array_equal(frcvr_target_D4, mg.at_node['flow__receiver_node'])
    assert_array_equal(upids_target_D4, mg.at_node['flow__upstream_node_order'])
    assert_array_equal(links2rcvr_target_D4,
                       mg.at_node['flow__link_to_receiver_node'])
    assert_array_almost_equal(steepest_target_D4,
                              mg.at_node['topographic__steepest_slope'])


@with_setup(setup_internal_closed)
def test_internal_closed():
    """Test closed nodes in the core of the grid."""
    fr = FlowRouter(mg)
    fr.route_flow()
    assert_array_almost_equal(A_target, mg.at_node['drainage_area'])
    assert_array_equal(frcvr_target, mg.at_node['flow__receiver_node'])
    assert_array_equal(links2rcvr_target, mg.at_node['flow__link_to_receiver_node'])
    assert_array_almost_equal(A_target, mg.at_node['water__discharge'])
    assert_array_almost_equal(steepest_target,
                              mg.at_node['topographic__steepest_slope'])


@with_setup(setup_voronoi)
def test_voronoi():
    """Test routing on a (radial) voronoi."""
    fr.route_flow()
    assert_array_almost_equal(vmg.at_node['drainage_area'][vmg.core_nodes],
                              A_target_core)
    assert_almost_equal(vmg.at_node['drainage_area'][12], A_target_outlet)


@with_setup(setup_voronoi_closedinternal)
def test_voronoi_closedinternal():
    """Test routing on a (radial) voronoi, but with a closed interior node."""
    fr.route_flow()
    # for i in range(vmg.number_of_nodes):
    #     print i, vmg.node_x[i], vmg.node_y[i], vmg.status_at_node[i], \
    #             vmg.at_node['drainage_area'][i], vmg.at_node['flow__receiver_node'][i], \
    #             vmg.at_node['topographic__elevation'][i]

    assert_array_almost_equal(vmg.at_node['drainage_area'][vmg.core_nodes],
                              A_target_internal)
    assert_almost_equal(vmg.at_node['drainage_area'][12], A_target_outlet)

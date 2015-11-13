# -*- coding: utf-8 -*-
"""
test_flow_routing:

Created on Thurs Nov 12, 2015

@author: dejh

This just tests the functionality of the router in toto - no attempt is made
to test the submodules.
"""

import landlab
from landlab import RasterModelGrid, FieldError
from landlab.components.flow_routing.route_flow_dn import FlowRouter
from numpy import sin, pi
import numpy as np  # for use of np.round
from numpy.testing import assert_array_equal, assert_array_almost_equal
from landlab import BAD_INDEX_VALUE as XX
from nose.tools import (with_setup, assert_true, assert_false, assert_raises,
                        assert_almost_equal, assert_equal)
try:
    from nose.tools import (assert_is, assert_set_equal, assert_dict_equal)
except ImportError:
    from landlab.testing.tools import (assert_is, assert_set_equal,
                                       assert_dict_equal)


def setup_dans_grid1():
    """
    Create a 5x5 test grid.
    This is a sheet flow test.
    """
    global fr, mg
    global (z, Q_in, A_target, frcvr_target, upids_target, Q_target,
            steepest_target, links2rcvr_target)

    mg = RasterModelGrid((5, 5), spacing=(10., 10.))

    z = mg.node_x.copy()

    Q_in = np.full(25, 2.)

    A_target = np.array([1.,  1.,  1.,  1.,  1.,
                         4.,  3.,  2.,  1.,  1.,
                         4.,  3.,  2.,  1.,  1.,
                         4.,  3.,  2.,  1.,  1.,
                         1.,  1.,  1.,  1.,  1.])*100.

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
    links2rcvr_target[mg.core_nodes] = np.array([24, 25, 26,
                                                 28, 29, 30,
                                                 32, 33, 34])

    Q_target = A_target * 2.  # only once Q_in is used

    steepest_target = np.array([0.,  0.,  0.,  0.,  0.,
                                0.,  1.,  1.,  1.,  0.,
                                0.,  1.,  1.,  1.,  0.,
                                0.,  1.,  1.,  1.,  0.,
                                0.,  0.,  0.,  0.,  0.])

    mg.add_field('node', 'topographic__elevation', z, units='-')


@with_setup(setup_dans_grid1)
def check_fields():
    """
    Check to make sure the right fields have been created.
    """
    fr = FlowRouter(mg)
    assert_array_equal(z, mg.at_node['topographic__elevation'])
    assert_array_equal(np.zeros(25), mg.at_node['drainage_area'])
    assert_array_equal(np.ones(25), mg.at_node['water__volume_flux_in'])

    fr = FlowRouter(mg, 'test_fr_input.txt')
    assert_array_equal(np.full(25, 2.), mg.at_node['water__volume_flux_in'])


@with_setup(setup_dans_grid1)
def check_field_input():
    """
    Check we can successfully pass water__volume_flux_in
    """
    mg.add_field('node', 'water__volume_flux_in',
                 np.full(25, 3.), units='m**3/s')
    fr = FlowRouter(mg)
    assert_array_equal(np.full(25, 3.), mg.at_node['water__volume_flux_in'])
    fr = FlowRouter(mg, 'test_fr_input.txt')
    assert_array_equal(np.full(25, 2.), mg.at_node['water__volume_flux_in'])


@with_setup(setup_dans_grid1)
def test_get_lake_ext_margin():
    lake = np.array([16, 17, 23, 24, 25, 30, 31, 32])
    ext_margin_returned = hf.get_lake_ext_margin(lake)
    ext_margin = np.array([8, 9, 10, 11, 15, 18, 19, 22, 26, 29, 33, 36, 37,
                           38, 39, 40])
    assert_array_equal(ext_margin_returned, ext_margin)


@with_setup(setup_dans_grid1)
def test_get_lake_int_margin():
    lake = np.array([16, 17, 18, 23, 24, 25, 26, 30, 31, 32])
    ext_margin = np.array([8, 9, 10, 11, 12, 15, 19, 20, 22, 27, 29, 33, 34,
                           36, 37, 38, 39, 40])
    int_margin_returned = hf.get_lake_int_margin(lake, ext_margin)
    int_margin = np.array([16, 17, 18, 23, 25, 26, 30, 31, 32])
    assert_array_equal(int_margin_returned, int_margin)


@with_setup(setup_dans_grid1)
def test_drainage_directions_change():
    lake = np.array([22, 23])
    old_elevs = np.ones(49, dtype=float)
    old_elevs[lake] = 0.
    new_elevs = old_elevs.copy()
    new_elevs[40] = 2.
    cond = hf.drainage_directions_change(lake, old_elevs, new_elevs)
    assert_false(cond)
    new_elevs[23] = 0.5
    cond = hf.drainage_directions_change(lake, old_elevs, new_elevs)
    assert_false(cond)
    new_elevs[23] = 1.
    cond = hf.drainage_directions_change(lake, old_elevs, new_elevs)
    assert_false(cond)
    new_elevs[23] = 1.2
    cond = hf.drainage_directions_change(lake, old_elevs, new_elevs)
    assert_true(cond)


@with_setup(setup_dans_grid1)
def test_add_slopes():
    new_z = z.copy()
    outlet_elev = z[outlet]
    hf._elev[lake] = outlet_elev
    rt2 = np.sqrt(2.)
    slope_to_add = 0.1
    lake_map = np.empty_like(z)
    lake_map.fill(XX)
    lake_map[lake] = lake_code
    hf._lf._lake_map = lake_map
    hf.lake_nodes_treated = np.array([], dtype=int)
    dists = mg.get_distances_of_nodes_to_point((mg.node_x[outlet],
                                                mg.node_y[outlet]))
    new_z[lake] = outlet_elev
    new_z[lake] += dists[lake]*slope_to_add
    # test the ones we can do easily analytically separately
    straight_north = np.array([23, 16])
    off_angle = 24
    elevs_out, lake_out = hf.add_slopes(slope_to_add, outlet, lake_code)
    assert_array_equal(slope_to_add*(np.arange(2.)+1.)+outlet_elev,
                       elevs_out[straight_north])
    assert_almost_equal(slope_to_add*rt2+outlet_elev, elevs_out[off_angle])
    assert_array_equal(new_z, elevs_out)
    assert_array_equal(lake, lake_out)


@with_setup(setup_dans_grid2)
def test_filler_flat():
    """
    Very simple, though possibly degerate, case, filling a 3x3 hole up to
    the flat surface surrounding it.
    """
    hf.fill_pits()
    assert_array_equal(hf._elev[lake], np.ones(9, dtype=float))
    assert_array_equal(mg.at_node['topographic__elevation'][lake],
                       np.ones(9, dtype=float))


@with_setup(setup_dans_grid3)
def test_filler_inclined():
    """
    Tests a flat fill into an inclined surface, with two holes.
    """
    hf.fill_pits()
    assert_array_equal(mg.at_node['topographic__elevation'][lake1],
                       np.ones(9, dtype=float)*4.)
    assert_array_equal(mg.at_node['topographic__elevation'][lake2],
                       np.ones(4, dtype=float)*7.)


@with_setup(setup_dans_grid3)
def test_filler_inclined2():
    """
    Tests an inclined fill into an inclined surface, with two holes.
    """
    z_init = z.copy()
    hf.fill_pits(apply_slope=True)
    hole1 = np.array([4.00009091, 4.00018182, 4.00027273, 4.00063636,
                      4.00045455, 4.00036364, 4.00081818, 4.00072727,
                      4.00054545])
    hole2 = np.array([7.16666667, 7.33333333, 7.66666667, 7.5])
    assert_array_almost_equal(mg.at_node['topographic__elevation'][lake1],
                              hole1)
    assert_array_almost_equal(mg.at_node['topographic__elevation'][lake2],
                              hole2)
    fr.route_flow()
    assert_equal(mg.at_node['flow_sinks'][mg.core_nodes].sum(), 0)


@with_setup(setup_dans_grid4)
def test_stupid_shaped_hole():
    """
    Tests inclined fill into a surface with a deliberately awkward shape.
    """
    hf.fill_pits(apply_slope=True)
    hole1 = np.array([4.00007692, 4.00015385, 4.00023077, 4.00053846,
                      4.00038462, 4.00030769, 4.00069231, 4.00061538,
                      4.00046154, 4.00076923, 4.00084615])
    hole2 = np.array([7.4, 7.2, 7.6])
    # print this to check out the funky drainage arrangement...
    # print(mg.at_node['topographic__elevation'].reshape((10, 10))[3:8, 4:7])
    assert_array_almost_equal(mg.at_node['topographic__elevation'][lake1],
                              hole1)
    assert_array_almost_equal(mg.at_node['topographic__elevation'][lake2],
                              hole2)
    fr.route_flow()
    assert_equal(mg.at_node['flow_sinks'][mg.core_nodes].sum(), 0)


@with_setup(setup_dans_grid5)
def test_D4_routing():
    """
    Tests inclined fill into a surface with a deliberately awkward shape.
    This is testing D4 routing.
    """
    hf.fill_pits(apply_slope=True)
    hole1 = np.array([4.00016667, 4.00025, 4.00033333, 4.00008333, 4.00041667,
                      4.0005, 4.00083333, 4.00066667, 4.00058333, 4.00075,
                      4.334])
    hole2 = np.array([7.6, 7.2, 7.4])
    # np.array([34, 35, 36, 44, 45, 46, 54, 55, 56, 65, 74])
    # print this to check out the funky drainage arrangement...
    # print(mg.at_node['topographic__elevation'].reshape((10, 10))[3:8, 4:7])
    assert_array_almost_equal(mg.at_node['topographic__elevation'][lake1],
                              hole1)
    assert_array_almost_equal(mg.at_node['topographic__elevation'][lake2],
                              hole2)
    fr.route_flow(method='D4')
    assert_equal(mg.at_node['flow_sinks'][mg.core_nodes].sum(), 0)


@with_setup(setup_dans_grid5)
def test_D4_filling():
    """
    Tests inclined fill into a surface with a deliberately awkward shape.
    This is testing D4 without inclining the surface.
    """
    hf.fill_pits(apply_slope=False)
    hole1 = 4.*np.ones_like(lake1, dtype=float)
    hole1[-1] += 0.001
    hole2 = 7.*np.ones_like(lake2, dtype=float)
    # np.array([34, 35, 36, 44, 45, 46, 54, 55, 56, 65, 74])
    # print this to check out the funky drainage arrangement...
    # print(mg.at_node['topographic__elevation'].reshape((10, 10))[3:8, 4:7])
    assert_array_almost_equal(mg.at_node['topographic__elevation'][lake1],
                              hole1)
    assert_array_almost_equal(mg.at_node['topographic__elevation'][lake2],
                              hole2)

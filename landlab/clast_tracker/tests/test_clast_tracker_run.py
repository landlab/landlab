# -*- coding: utf-8 -*-
"""
Unit tests for landlab.clast_tracker
Tests run_on_step output

Last updated 02/11/2019

"""

import pytest
import numpy as np
from landlab.components import FlowAccumulator #, LinearDiffuser



def test_run_one_step_flat(cc_flat, grid_flat):
    fa = FlowAccumulator(grid_flat,
                    'topographic__elevation',
                    flow_director='FlowDirectorSteepest')
    fa.run_one_step()
    kappa = 0.001
    cc_flat.run_one_step(dt=10.,
                         Si=1.2,
                         kappa=kappa,
                         uplift=None,
                         hillslope_river__threshold=1e4,
                         lateral_spreading='off',
                         disturbance_fqcy=1.,
                         d_star=1.)

    fa.run_one_step()
    cc_flat.run_one_step(dt=10.,
                          Si=1.2,
                          kappa=kappa,
                          uplift=None,
                          hillslope_river__threshold=1e4,
                          lateral_spreading='off',
                          disturbance_fqcy=1.,
                          d_star=1.)

    assert cc_flat.slope__WE.values[0] == 0.
    assert np.isnan(cc_flat.slope__WE.values[1])
    assert cc_flat.slope__SN.values[0] == 0.
    assert np.isnan(cc_flat.slope__SN.values[1])
    assert np.all(np.isnan(cc_flat.slope__steepest_azimuth.values[0]))
    assert np.allclose(cc_flat.slope__steepest_dip.values, [0., 0.])
    assert np.allclose(cc_flat.target_node_flag.values, [-1, -1])
    assert np.allclose(cc_flat.target_node.values, [12, 13])
    assert cc_flat.clast__x.values[0,1] == cc_flat.clast__x.values[0,0]
    assert cc_flat.clast__y.values[0,1] == cc_flat.clast__y.values[0,0]
    assert cc_flat.clast__elev.values[0,1] == cc_flat.clast__elev.values[0,0]
    assert cc_flat.clast__elev.values[1,1] == cc_flat.clast__elev.values[1,0]
    assert np.allclose(cc_flat.total_travelled_dist.values, [0., 0.])

def test_run_one_step_south(cc_south, grid_south):
    fa = FlowAccumulator(grid_south,
                         'topographic__elevation',
                         flow_director='FlowDirectorSteepest')
    kappa = 0.001
#not used
#    ld = LinearDiffuser(grid_south,
#                        linear_diffusivity=kappa)

    fa.run_one_step()
    cc_south.run_one_step(dt=10.,
                          Si=1.2,
                          kappa=kappa,
                          uplift=None,
                          hillslope_river__threshold=1e4,
                          lateral_spreading='off',
                          disturbance_fqcy=1.,
                          d_star=1.)

    # attributes
    assert cc_south.attrs['kappa'] == 0.001
    # coordinates
    assert cc_south.time_coordinates == [0.0, 10.0]
    # slope
    assert np.isclose(cc_south.slope__steepest_azimuth.values[0], 4.71238898)
    assert np.isclose(cc_south.slope__steepest_dip.values[0], 0.2914567944)
    # travel distance
    assert (cc_south.total_travelled_dist.values[0] ==
            cc_south.hop_length.values[0,1])

def test_run_one_step_east(cc_east, grid_east):
    fa = FlowAccumulator(grid_east,
                         'topographic__elevation',
                         flow_director='FlowDirectorSteepest')
    kappa = 0.001

    fa.run_one_step()
    cc_east.run_one_step(dt=10.,
                          Si=1.2,
                          kappa=kappa,
                          uplift=None,
                          hillslope_river__threshold=1e4,
                          lateral_spreading='off',
                          disturbance_fqcy=1.,
                          d_star=1.)

    assert np.allclose(cc_east.slope__steepest_azimuth.values, [0.0, 0.0])
    assert np.allclose(cc_east.slope__steepest_dip.values,
                       [0.2914567944, 0.2914567944])
    assert np.allclose(cc_east.total_travelled_dist.values,
                       cc_east.hop_length.values[:,1])
    assert np.all(cc_east.hop_length.values[:,1] > 0)
    assert cc_east.close2boundary.values[1] == 1.0

def test_run_one_step_north(cc_north, grid_north):
    fa = FlowAccumulator(grid_north,
                         'topographic__elevation',
                         flow_director='FlowDirectorSteepest')
    kappa = 0.001

    fa.run_one_step()
    cc_north.run_one_step(dt=10.,
                          Si=1.2,
                          kappa=kappa,
                          uplift=None,
                          hillslope_river__threshold=1e4,
                          lateral_spreading='off',
                          disturbance_fqcy=1.,
                          d_star=1.)

    assert np.allclose(cc_north.slope__steepest_azimuth.values,
                       [1.57079633, 1.57079633])
    assert np.allclose(cc_north.slope__steepest_dip.values,
                       [0.2914567944, 0.2914567944])
    assert np.allclose(cc_north.total_travelled_dist.values,
                       cc_north.hop_length.values[:,1])
    assert np.all(cc_north.hop_length.values[:,1] > 0)
    assert cc_north.close2boundary.values[1] == 1.0
    assert np.isnan(cc_north.close2boundary.values[0])

def test_run_one_step_west(cc_west, grid_west):
    fa = FlowAccumulator(grid_west,
                         'topographic__elevation',
                         flow_director='FlowDirectorSteepest')
    kappa = 0.001

    fa.run_one_step()
    cc_west.run_one_step(dt=10.,
                          Si=1.2,
                          kappa=kappa,
                          uplift=None,
                          hillslope_river__threshold=1e4,
                          lateral_spreading='off',
                          disturbance_fqcy=1.,
                          d_star=1.)

    assert np.allclose(cc_west.slope__steepest_azimuth.values,
                       [3.1415926535, 3.1415926535])
    assert np.allclose(cc_west.slope__steepest_dip.values,
                       [0.2914567944, 0.2914567944])
    assert np.allclose(cc_west.total_travelled_dist.values,
                       cc_west.hop_length.values[:,1])
    assert np.all(cc_west.hop_length.values[:,1] > 0)
    assert cc_west.close2boundary.values[1] == 1.0
    assert np.isnan(cc_west.close2boundary.values[0])

# Diagonals

def test_run_one_step_ne(cc_ne, grid_ne):
    fa = FlowAccumulator(grid_ne,
                         'topographic__elevation',
                         flow_director='D8')
    kappa = 0.001

    fa.run_one_step()
    cc_ne.run_one_step(dt=10.,
                          Si=1.2,
                          kappa=kappa,
                          uplift=None,
                          hillslope_river__threshold=1e4,
                          lateral_spreading='off',
                          disturbance_fqcy=1.,
                          d_star=1.)

    assert np.allclose(cc_ne.target_node.values, [17, 18, 18])
    assert np.allclose(cc_ne.slope__steepest_azimuth.values,
                       [0.785398, 0.785398, 1.570796])
    assert np.allclose(cc_ne.slope__steepest_dip.values,
                       [0.401247, 0.401247, 0.291457])
    assert np.allclose(cc_ne.total_travelled_dist.values,
                       cc_ne.hop_length.values[:,1])
    assert np.all(cc_ne.hop_length.values[:,1] > 0)
    assert cc_ne.close2boundary.values[0] == 1.0
    assert np.isnan(cc_ne.close2boundary.values[1])

def test_run_one_step_nw(cc_nw, grid_nw):
    fa = FlowAccumulator(grid_nw,
                         'topographic__elevation',
                         flow_director='D8')
    kappa = 0.001

    fa.run_one_step()
    cc_nw.run_one_step(dt=10.,
                          Si=1.2,
                          kappa=kappa,
                          uplift=None,
                          hillslope_river__threshold=1e4,
                          lateral_spreading='off',
                          disturbance_fqcy=1.,
                          d_star=1.)

    assert np.allclose(cc_nw.target_node.values, [16, 16, 17])
    assert np.allclose(cc_nw.slope__steepest_azimuth.values,
                       [1.570796, 2.356194, 2.356194])
    assert np.allclose(cc_nw.slope__steepest_dip.values,
                       [0.291457, 0.401247, 0.401247])
    assert np.allclose(cc_nw.total_travelled_dist.values,
                       cc_nw.hop_length.values[:,1])
    assert np.all(cc_nw.hop_length.values[:,1] > 0)
    assert cc_nw.close2boundary.values[0] == 1.0
    assert np.isnan(cc_nw.close2boundary.values[1])

def test_run_one_step_sw(cc_sw, grid_sw):
    fa = FlowAccumulator(grid_sw,
                         'topographic__elevation',
                         flow_director='D8')
    kappa = 0.001

    fa.run_one_step()
    cc_sw.run_one_step(dt=10.,
                          Si=1.2,
                          kappa=kappa,
                          uplift=None,
                          hillslope_river__threshold=1e4,
                          lateral_spreading='off',
                          disturbance_fqcy=1.,
                          d_star=1.)

    assert np.allclose(cc_sw.target_node.values, [6, 6, 7])
    assert np.allclose(cc_sw.slope__steepest_azimuth.values,
                       [4.712389, 3.926991, 3.926991])
    assert np.allclose(cc_sw.slope__steepest_dip.values,
                       [0.291457, 0.401247, 0.401247])
    assert np.allclose(cc_sw.total_travelled_dist.values,
                       cc_sw.hop_length.values[:,1])
    assert np.all(cc_sw.hop_length.values[:,1] > 0)
    assert cc_sw.close2boundary.values[0] == 1.0
    assert np.isnan(cc_sw.close2boundary.values[1])

def test_run_one_step_se(cc_se, grid_se):
    fa = FlowAccumulator(grid_se,
                         'topographic__elevation',
                         flow_director='D8')
    kappa = 0.001

    fa.run_one_step()
    cc_se.run_one_step(dt=10.,
                          Si=1.2,
                          kappa=kappa,
                          uplift=None,
                          hillslope_river__threshold=1e4,
                          lateral_spreading='off',
                          disturbance_fqcy=1.,
                          d_star=1.)

    assert np.allclose(cc_se.target_node.values, [7, 8, 8])
    assert np.allclose(cc_se.slope__steepest_azimuth.values,
                       [5.497787, 5.497787, 4.712389])
    assert np.allclose(cc_se.slope__steepest_dip.values,
                       [0.401247, 0.401247, 0.291457])
    assert np.allclose(cc_se.total_travelled_dist.values,
                       cc_se.hop_length.values[:,1])
    assert np.all(cc_se.hop_length.values[:,1] > 0)
    assert cc_se.close2boundary.values[0] == 1.0
    assert np.isnan(cc_se.close2boundary.values[1])



def test_run_one_with_dev(cc_south, grid_south):
    fa = FlowAccumulator(grid_south,
                         'topographic__elevation',
                         flow_director='FlowDirectorSteepest')
    kappa = 0.001
    fa.run_one_step()
    cc_south.run_one_step(dt=10.,
                                      Si=1.2,
                                      kappa=kappa,
                                      uplift=None,
                                      hillslope_river__threshold=1e4,
                                      lateral_spreading='on',
                                      disturbance_fqcy=1.,
                                      d_star=1.)
    assert cc_south.change_x.values[0] != 0
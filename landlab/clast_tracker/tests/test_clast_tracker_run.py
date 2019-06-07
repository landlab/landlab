# -*- coding: utf-8 -*-
"""
Unit tests for landlab.clast_tracker
Tests run_on_step output

Last updated 02/11/2019

"""

import pytest
import numpy as np
from landlab.components import FlowAccumulator

kappa = 0.001


def test_run_one_step_flat(cc_flat, grid_flat):
    fa = FlowAccumulator(
        grid_flat,
        "topographic__elevation",
        flow_director="FlowDirectorSteepest",
    )
    fa.run_one_step()
    cc_flat.run_one_step(
        dt=10.0,
        Si=1.2,
        kappa=kappa,
        uplift=None,
        hillslope_river__threshold=1e4,
        lateral_spreading="off",
        disturbance_fqcy=1.0,
        d_star=1.0,
    )

    fa.run_one_step()
    cc_flat.run_one_step(
        dt=10.0,
        Si=1.2,
        kappa=kappa,
        uplift=None,
        hillslope_river__threshold=1e4,
        lateral_spreading="off",
        disturbance_fqcy=1.0,
        d_star=1.0,
    )

    assert cc_flat.slope__WE.values[0] == 0.0
    assert np.isnan(cc_flat.slope__WE.values[1])
    assert cc_flat.slope__SN.values[0] == 0.0
    assert np.isnan(cc_flat.slope__SN.values[1])
    assert np.all(
        np.isnan(cc_flat.slope__steepest_azimuth.values[0])
    )
    assert np.allclose(
        cc_flat.slope__steepest_dip.values, [0.0, 0.0]
    )
    assert np.allclose(
        cc_flat.target_node_flag.values, [-1, -1]
    )
    assert np.allclose(cc_flat.target_node.values, [12, 13])
    assert (
        cc_flat.clast__x.values[0, 1]
        == cc_flat.clast__x.values[0, 0]
    )
    assert (
        cc_flat.clast__y.values[0, 1]
        == cc_flat.clast__y.values[0, 0]
    )
    assert (
        cc_flat.clast__elev.values[0, 1]
        == cc_flat.clast__elev.values[0, 0]
    )
    assert (
        cc_flat.clast__elev.values[1, 1]
        == cc_flat.clast__elev.values[1, 0]
    )
    assert np.allclose(
        cc_flat.total_travelled_dist.values, [0.0, 0.0]
    )


def test_run_one_step_south(cc_south, grid_south):
    fa = FlowAccumulator(
        grid_south,
        "topographic__elevation",
        flow_director="FlowDirectorSteepest",
    )
    fa.run_one_step()
    cc_south.run_one_step(
        dt=10.0,
        Si=1.2,
        kappa=kappa,
        uplift=None,
        hillslope_river__threshold=1e4,
        lateral_spreading="off",
        disturbance_fqcy=1.0,
        d_star=1.0,
    )
    # attributes
    assert cc_south.attrs["kappa"] == 0.001
    # coordinates
    assert cc_south.time_coordinates == [0.0, 10.0]
    # slope
    assert np.isclose(
        cc_south.slope__steepest_azimuth.values[0],
        4.71238898,
    )
    assert np.isclose(
        cc_south.slope__steepest_dip.values[0], 0.2914567944
    )
    # travel distance
    assert (
        cc_south.total_travelled_dist.values[0]
        == cc_south.hop_length.values[0, 1]
    )


def test_run_one_step_east(cc_east, grid_east):
    fa = FlowAccumulator(
        grid_east,
        "topographic__elevation",
        flow_director="FlowDirectorSteepest",
    )
    fa.run_one_step()
    cc_east.run_one_step(
        dt=10.0,
        Si=1.2,
        kappa=kappa,
        uplift=None,
        hillslope_river__threshold=1e4,
        lateral_spreading="off",
        disturbance_fqcy=1.0,
        d_star=1.0,
    )

    assert np.allclose(
        cc_east.slope__steepest_azimuth.values, [0.0, 0.0]
    )
    assert np.allclose(
        cc_east.slope__steepest_dip.values,
        [0.2914567944, 0.2914567944],
    )
    assert np.allclose(
        cc_east.total_travelled_dist.values,
        cc_east.hop_length.values[:, 1],
    )
    assert np.all(cc_east.hop_length.values[:, 1] > 0)
    assert cc_east.close2boundary.values[1] == 1.0


def test_run_one_step_north(cc_north, grid_north):
    fa = FlowAccumulator(
        grid_north,
        "topographic__elevation",
        flow_director="FlowDirectorSteepest",
    )

    fa.run_one_step()
    cc_north.run_one_step(
        dt=10.0,
        Si=1.2,
        kappa=kappa,
        uplift=None,
        hillslope_river__threshold=1e4,
        lateral_spreading="off",
        disturbance_fqcy=1.0,
        d_star=1.0,
    )

    assert np.allclose(
        cc_north.slope__steepest_azimuth.values,
        [1.57079633, 1.57079633],
    )
    assert np.allclose(
        cc_north.slope__steepest_dip.values,
        [0.2914567944, 0.2914567944],
    )
    assert np.allclose(
        cc_north.total_travelled_dist.values,
        cc_north.hop_length.values[:, 1],
    )
    assert np.all(cc_north.hop_length.values[:, 1] > 0)
    assert cc_north.close2boundary.values[1] == 1.0
    assert np.isnan(cc_north.close2boundary.values[0])


def test_run_one_step_west(cc_west, grid_west):
    fa = FlowAccumulator(
        grid_west,
        "topographic__elevation",
        flow_director="FlowDirectorSteepest",
    )
    fa.run_one_step()
    cc_west.run_one_step(
        dt=10.0,
        Si=1.2,
        kappa=kappa,
        uplift=None,
        hillslope_river__threshold=1e4,
        lateral_spreading="off",
        disturbance_fqcy=1.0,
        d_star=1.0,
    )
    assert np.allclose(
        cc_west.slope__steepest_azimuth.values,
        [3.1415926535, 3.1415926535],
    )
    assert np.allclose(
        cc_west.slope__steepest_dip.values,
        [0.2914567944, 0.2914567944],
    )
    assert np.allclose(
        cc_west.total_travelled_dist.values,
        cc_west.hop_length.values[:, 1],
    )
    assert np.all(cc_west.hop_length.values[:, 1] > 0)
    assert cc_west.close2boundary.values[1] == 1.0
    assert np.isnan(cc_west.close2boundary.values[0])


# Diagonals
def test_run_one_step_ne(cc_ne, grid_ne):
    fa = FlowAccumulator(
        grid_ne,
        "topographic__elevation",
        flow_director="D8",
    )
    fa.run_one_step()
    cc_ne.run_one_step(
        dt=10.0,
        Si=1.2,
        kappa=kappa,
        uplift=None,
        hillslope_river__threshold=1e4,
        lateral_spreading="off",
        disturbance_fqcy=1.0,
        d_star=1.0,
    )
    assert np.allclose(
        cc_ne.target_node.values, [17, 18, 18]
    )
    assert np.allclose(
        cc_ne.slope__steepest_azimuth.values,
        [0.785398, 0.785398, 1.570796],
    )
    assert np.allclose(
        cc_ne.slope__steepest_dip.values,
        [0.401247, 0.401247, 0.291457],
    )
    assert np.allclose(
        cc_ne.total_travelled_dist.values,
        cc_ne.hop_length.values[:, 1],
    )
    assert np.all(cc_ne.hop_length.values[:, 1] > 0)
    assert cc_ne.close2boundary.values[0] == 1.0
    assert np.isnan(cc_ne.close2boundary.values[1])


def test_run_one_step_nw(cc_nw, grid_nw):
    fa = FlowAccumulator(
        grid_nw,
        "topographic__elevation",
        flow_director="D8",
    )
    fa.run_one_step()
    cc_nw.run_one_step(
        dt=10.0,
        Si=1.2,
        kappa=kappa,
        uplift=None,
        hillslope_river__threshold=1e4,
        lateral_spreading="off",
        disturbance_fqcy=1.0,
        d_star=1.0,
    )
    assert np.allclose(
        cc_nw.target_node.values, [16, 16, 17]
    )
    assert np.allclose(
        cc_nw.slope__steepest_azimuth.values,
        [1.570796, 2.356194, 2.356194],
    )
    assert np.allclose(
        cc_nw.slope__steepest_dip.values,
        [0.291457, 0.401247, 0.401247],
    )
    assert np.allclose(
        cc_nw.total_travelled_dist.values,
        cc_nw.hop_length.values[:, 1],
    )
    assert np.all(cc_nw.hop_length.values[:, 1] > 0)
    assert cc_nw.close2boundary.values[0] == 1.0
    assert np.isnan(cc_nw.close2boundary.values[1])


def test_run_one_step_sw(cc_sw, grid_sw):
    fa = FlowAccumulator(
        grid_sw,
        "topographic__elevation",
        flow_director="D8",
    )
    fa.run_one_step()
    cc_sw.run_one_step(
        dt=10.0,
        Si=1.2,
        kappa=kappa,
        uplift=None,
        hillslope_river__threshold=1e4,
        lateral_spreading="off",
        disturbance_fqcy=1.0,
        d_star=1.0,
    )
    assert np.allclose(cc_sw.target_node.values, [6, 6, 7])
    assert np.allclose(
        cc_sw.slope__steepest_azimuth.values,
        [4.712389, 3.926991, 3.926991],
    )
    assert np.allclose(
        cc_sw.slope__steepest_dip.values,
        [0.291457, 0.401247, 0.401247],
    )
    assert np.allclose(
        cc_sw.total_travelled_dist.values,
        cc_sw.hop_length.values[:, 1],
    )
    assert np.all(cc_sw.hop_length.values[:, 1] > 0)
    assert cc_sw.close2boundary.values[0] == 1.0
    assert np.isnan(cc_sw.close2boundary.values[1])


def test_run_one_step_se(cc_se, grid_se):
    fa = FlowAccumulator(
        grid_se,
        "topographic__elevation",
        flow_director="D8",
    )
    fa.run_one_step()
    cc_se.run_one_step(
        dt=10.0,
        Si=1.2,
        kappa=kappa,
        uplift=None,
        hillslope_river__threshold=1e4,
        lateral_spreading="off",
        disturbance_fqcy=1.0,
        d_star=1.0,
    )
    assert np.allclose(cc_se.target_node.values, [7, 8, 8])
    assert np.allclose(
        cc_se.slope__steepest_azimuth.values,
        [5.497787, 5.497787, 4.712389],
    )
    assert np.allclose(
        cc_se.slope__steepest_dip.values,
        [0.401247, 0.401247, 0.291457],
    )
    assert np.allclose(
        cc_se.total_travelled_dist.values,
        cc_se.hop_length.values[:, 1],
    )
    assert np.all(cc_se.hop_length.values[:, 1] > 0)
    assert cc_se.close2boundary.values[0] == 1.0
    assert np.isnan(cc_se.close2boundary.values[1])


# Intermediate slopes


def test_run_one_step_ene(cc_ene, grid_ene):
    fa = FlowAccumulator(
        grid_ene,
        "topographic__elevation",
        flow_director="D8",
    )
    fa.run_one_step()
    cc_ene.run_one_step(
        dt=10.0,
        Si=1.2,
        kappa=kappa,
        uplift=None,
        hillslope_river__threshold=1e4,
        lateral_spreading="off",
        disturbance_fqcy=1.0,
        d_star=1.0,
    )
    assert np.allclose(
        cc_ene.target_node.values, [17, 13, 18]
    )
    assert np.allclose(
        cc_ene.slope__steepest_azimuth.values,
        [0.78539816, 0.46364761, 1.57079633],
    )
    assert np.allclose(
        cc_ene.slope__steepest_dip.values,
        [0.56675233, 0.59087275, 0.29145679],
    )
    assert np.allclose(
        cc_ene.total_travelled_dist.values,
        cc_ene.hop_length.values[:, 1],
    )
    assert np.all(cc_ene.hop_length.values[:, 1] > 0)
    assert cc_ene.close2boundary.values[0] == 1.0
    assert np.isnan(cc_ene.close2boundary.values[1])


def test_run_one_step_nne(cc_nne, grid_nne):
    fa = FlowAccumulator(
        grid_nne,
        "topographic__elevation",
        flow_director="D8",
    )
    fa.run_one_step()
    cc_nne.run_one_step(
        dt=10.0,
        Si=1.2,
        kappa=kappa,
        uplift=None,
        hillslope_river__threshold=1e4,
        lateral_spreading="off",
        disturbance_fqcy=1.0,
        d_star=1.0,
    )
    assert np.allclose(
        cc_nne.target_node.values, [17, 17, 18]
    )
    assert np.allclose(
        cc_nne.slope__steepest_azimuth.values,
        [0.78539816, 1.10714872, 1.57079633],
    )
    assert np.allclose(
        cc_nne.slope__steepest_dip.values,
        [0.56675233, 0.59087275, 0.5404195],
    )
    assert np.allclose(
        cc_nne.total_travelled_dist.values,
        cc_nne.hop_length.values[:, 1],
    )
    assert np.all(cc_nne.hop_length.values[:, 1] > 0)
    assert cc_nne.close2boundary.values[0] == 1.0
    assert np.isnan(cc_nne.close2boundary.values[1])


def test_run_one_step_nnw(cc_nnw, grid_nnw):
    fa = FlowAccumulator(
        grid_nnw,
        "topographic__elevation",
        flow_director="D8",
    )
    fa.run_one_step()
    cc_nnw.run_one_step(
        dt=10.0,
        Si=1.2,
        kappa=kappa,
        uplift=None,
        hillslope_river__threshold=1e4,
        lateral_spreading="off",
        disturbance_fqcy=1.0,
        d_star=1.0,
    )
    assert np.allclose(
        cc_nnw.target_node.values, [16, 17, 17]
    )
    assert np.allclose(
        cc_nnw.slope__steepest_azimuth.values,
        [1.57079633, 2.03444394, 2.35619449],
    )
    assert np.allclose(
        cc_nnw.slope__steepest_dip.values,
        [0.5404195, 0.59087275, 0.56675233],
    )
    assert np.allclose(
        cc_nnw.total_travelled_dist.values,
        cc_nnw.hop_length.values[:, 1],
    )
    assert np.all(cc_nnw.hop_length.values[:, 1] > 0)
    assert cc_nnw.close2boundary.values[0] == 1.0
    assert np.isnan(cc_nnw.close2boundary.values[1])


def test_run_one_step_wnw(cc_wnw, grid_wnw):
    fa = FlowAccumulator(
        grid_wnw,
        "topographic__elevation",
        flow_director="D8",
    )
    fa.run_one_step()
    cc_wnw.run_one_step(
        dt=10.0,
        Si=1.2,
        kappa=kappa,
        uplift=None,
        hillslope_river__threshold=1e4,
        lateral_spreading="off",
        disturbance_fqcy=1.0,
        d_star=1.0,
    )
    assert np.allclose(
        cc_wnw.target_node.values, [16, 11, 17]
    )
    assert np.allclose(
        cc_wnw.slope__steepest_azimuth.values,
        [1.57079633, 2.67794504, 2.35619449],
    )
    assert np.allclose(
        cc_wnw.slope__steepest_dip.values,
        [0.29145679, 0.59087275, 0.56675233],
    )
    assert np.allclose(
        cc_wnw.total_travelled_dist.values,
        cc_wnw.hop_length.values[:, 1],
    )
    assert np.all(cc_wnw.hop_length.values[:, 1] > 0)
    assert cc_wnw.close2boundary.values[0] == 1.0
    assert np.isnan(cc_wnw.close2boundary.values[1])


def test_run_one_step_wsw(cc_wsw, grid_wsw):
    fa = FlowAccumulator(
        grid_wsw,
        "topographic__elevation",
        flow_director="D8",
    )
    fa.run_one_step()
    cc_wsw.run_one_step(
        dt=10.0,
        Si=1.2,
        kappa=kappa,
        uplift=None,
        hillslope_river__threshold=1e4,
        lateral_spreading="off",
        disturbance_fqcy=1.0,
        d_star=1.0,
    )
    assert np.allclose(
        cc_wsw.target_node.values, [6, 11, 7]
    )
    assert np.allclose(
        cc_wsw.slope__steepest_azimuth.values,
        [4.71238898, 3.60524026, 3.92699082],
    )
    assert np.allclose(
        cc_wsw.slope__steepest_dip.values,
        [0.29145679, 0.59087275, 0.56675233],
    )
    assert np.allclose(
        cc_wsw.total_travelled_dist.values,
        cc_wsw.hop_length.values[:, 1],
    )
    assert np.all(cc_wsw.hop_length.values[:, 1] > 0)
    assert cc_wsw.close2boundary.values[0] == 1.0
    assert np.isnan(cc_wsw.close2boundary.values[1])


def test_run_one_step_ssw(cc_ssw, grid_ssw):
    fa = FlowAccumulator(
        grid_ssw,
        "topographic__elevation",
        flow_director="D8",
    )
    fa.run_one_step()
    cc_ssw.run_one_step(
        dt=10.0,
        Si=1.2,
        kappa=kappa,
        uplift=None,
        hillslope_river__threshold=1e4,
        lateral_spreading="off",
        disturbance_fqcy=1.0,
        d_star=1.0,
    )
    assert np.allclose(cc_ssw.target_node.values, [6, 7, 7])
    assert np.allclose(
        cc_ssw.slope__steepest_azimuth.values,
        [4.71238898, 4.24874137, 3.92699082],
    )
    assert np.allclose(
        cc_ssw.slope__steepest_dip.values,
        [0.5404195, 0.59087275, 0.56675233],
    )
    assert np.allclose(
        cc_ssw.total_travelled_dist.values,
        cc_ssw.hop_length.values[:, 1],
    )
    assert np.all(cc_ssw.hop_length.values[:, 1] > 0)
    assert cc_ssw.close2boundary.values[0] == 1.0
    assert np.isnan(cc_ssw.close2boundary.values[1])


def test_run_one_step_sse(cc_sse, grid_sse):
    fa = FlowAccumulator(
        grid_sse,
        "topographic__elevation",
        flow_director="D8",
    )
    fa.run_one_step()
    cc_sse.run_one_step(
        dt=10.0,
        Si=1.2,
        kappa=kappa,
        uplift=None,
        hillslope_river__threshold=1e4,
        lateral_spreading="off",
        disturbance_fqcy=1.0,
        d_star=1.0,
    )
    assert np.allclose(cc_sse.target_node.values, [7, 7, 8])
    assert np.allclose(
        cc_sse.slope__steepest_azimuth.values,
        [5.49778714, 5.17603659, 4.71238898],
    )
    assert np.allclose(
        cc_sse.slope__steepest_dip.values,
        [0.56675233, 0.59087275, 0.5404195],
    )
    assert np.allclose(
        cc_sse.total_travelled_dist.values,
        cc_sse.hop_length.values[:, 1],
    )
    assert np.all(cc_sse.hop_length.values[:, 1] > 0)
    assert cc_sse.close2boundary.values[0] == 1.0
    assert np.isnan(cc_sse.close2boundary.values[1])


def test_run_one_step_ese(cc_ese, grid_ese):
    fa = FlowAccumulator(
        grid_ese,
        "topographic__elevation",
        flow_director="D8",
    )
    fa.run_one_step()
    cc_ese.run_one_step(
        dt=10.0,
        Si=1.2,
        kappa=kappa,
        uplift=None,
        hillslope_river__threshold=1e4,
        lateral_spreading="off",
        disturbance_fqcy=1.0,
        d_star=1.0,
    )
    assert np.allclose(
        cc_ese.target_node.values, [7, 13, 8]
    )
    assert np.allclose(
        cc_ese.slope__steepest_azimuth.values,
        [5.49778714, 5.8195377, 4.71238898],
    )
    assert np.allclose(
        cc_ese.slope__steepest_dip.values,
        [0.56675233, 0.59087275, 0.29145679],
    )
    assert np.allclose(
        cc_ese.total_travelled_dist.values,
        cc_ese.hop_length.values[:, 1],
    )
    assert np.all(cc_ese.hop_length.values[:, 1] > 0)
    assert cc_ese.close2boundary.values[0] == 1.0
    assert np.isnan(cc_ese.close2boundary.values[1])


def test_run_one_with_dev(cc_south, grid_south):
    fa = FlowAccumulator(
        grid_south,
        "topographic__elevation",
        flow_director="FlowDirectorSteepest",
    )
    fa.run_one_step()
    cc_south.run_one_step(
        dt=10.0,
        Si=1.2,
        kappa=kappa,
        uplift=None,
        hillslope_river__threshold=1e4,
        lateral_spreading="on",
        disturbance_fqcy=1.0,
        d_star=1.0,
    )
    assert cc_south.change_x.values[0] != 0


def test_run_one_step_south_lowSc(cc_south, grid_south):
    fa = FlowAccumulator(
        grid_south,
        "topographic__elevation",
        flow_director="FlowDirectorSteepest",
    )
    fa.run_one_step()
    cc_south.run_one_step(
        dt=10.0,
        Si=0.1,
        kappa=kappa,
        uplift=None,
        hillslope_river__threshold=1e4,
        lateral_spreading="off",
        disturbance_fqcy=1.0,
        d_star=1.0,
    )
    assert cc_south.phantom(0)
    assert np.allclose(cc_south.clast__node.values, [2, 3])
    assert np.allclose(
        cc_south.total_travelled_dist,
        [1.56604598, 1.56604598],
    )


def test_run_one_step_south_river_highSc(
    cc_south, grid_south
):
    fa = FlowAccumulator(
        grid_south,
        "topographic__elevation",
        flow_director="FlowDirectorSteepest",
    )
    fa.run_one_step()
    cc_south.run_one_step(
        dt=10.0,
        Si=0.3,
        kappa=kappa,
        uplift=None,
        hillslope_river__threshold=1,
        lateral_spreading="off",
        disturbance_fqcy=1.0,
        d_star=1.0,
    )

    assert cc_south._cell_is_hillslope(0) is False
    assert np.allclose(cc_south.clast__node.values, [7, 8])
    assert np.allclose(
        cc_south.total_travelled_dist.values,
        [0.522015, 0.522015],
    )


def test_run_one_step_south_river_lowSc(
    cc_south, grid_south
):
    fa = FlowAccumulator(
        grid_south,
        "topographic__elevation",
        flow_director="FlowDirectorSteepest",
    )
    fa.run_one_step()
    cc_south.run_one_step(
        dt=10.0,
        Si=0.2,
        kappa=kappa,
        uplift=None,
        hillslope_river__threshold=1,
        lateral_spreading="off",
        disturbance_fqcy=1.0,
        d_star=1.0,
    )

    assert cc_south._cell_is_hillslope(0) is False
    assert np.allclose(cc_south.clast__node.values, [2, 3])
    assert np.allclose(
        cc_south.total_travelled_dist.values,
        [1.56604598, 1.56604598],
    )


def test_run_one_south_radius0(cc_south, grid_south):
    fa = FlowAccumulator(
        grid_south,
        "topographic__elevation",
        flow_director="FlowDirectorSteepest",
    )
    fa.run_one_step()
    cc_south.clast__radius.values[1, 0] = 0.0
    cc_south.run_one_step(
        dt=10.0,
        Si=1.2,
        kappa=kappa,
        uplift=None,
        hillslope_river__threshold=1e4,
        lateral_spreading="on",
        disturbance_fqcy=1.0,
        d_star=1.0,
    )

    assert cc_south.phantom(1)
    assert np.allclose(
        cc_south.total_travelled_dist.values,
        [0.00239766, 0.0],
    )


def test_run_one_uplift(cc_south, grid_south):
    fa = FlowAccumulator(
        grid_south,
        "topographic__elevation",
        flow_director="FlowDirectorSteepest",
    )
    fa.run_one_step()
    grid_south.add_field(
        "node",
        "spatially_var_uplift_field",
        [
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.1,
            0.2,
            0.3,
            0.0,
            0.0,
            0.1,
            0.2,
            0.3,
            0.0,
            0.0,
            0.1,
            0.2,
            0.3,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
        ],
    )
    uplift = "spatially_var_uplift_field"
    cc_south.run_one_step(
        dt=10.0,
        Si=1.2,
        kappa=kappa,
        uplift=uplift,
        hillslope_river__threshold=1e4,
        lateral_spreading="on",
        disturbance_fqcy=1.0,
        d_star=1.0,
    )
    grid_south.at_node["topographic__elevation"] += (
        10
        * grid_south.at_node["spatially_var_uplift_field"]
    )

    assert np.allclose(
        cc_south.uplift,
        [
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.1,
            0.2,
            0.3,
            0.0,
            0.0,
            0.1,
            0.2,
            0.3,
            0.0,
            0.0,
            0.1,
            0.2,
            0.3,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
        ],
    )
    assert np.allclose(
        cc_south.clast__elev.values[:, -1],
        grid_south.at_node["topographic__elevation"][
            cc_south.clast__node.values
        ],
    )

    uplift = 0.02
    cc_south.run_one_step(
        dt=10.0,
        Si=1.2,
        kappa=kappa,
        uplift=uplift,
        hillslope_river__threshold=1e4,
        lateral_spreading="on",
        disturbance_fqcy=1.0,
        d_star=1.0,
    )
    grid_south.at_node["topographic__elevation"] += (
        10 * uplift
    )

    assert np.allclose(
        cc_south.uplift,
        [
            0.02,
            0.02,
            0.02,
            0.02,
            0.02,
            0.02,
            0.02,
            0.02,
            0.02,
            0.02,
            0.02,
            0.02,
            0.02,
            0.02,
            0.02,
            0.02,
            0.02,
            0.02,
            0.02,
            0.02,
            0.02,
            0.02,
            0.02,
            0.02,
            0.02,
        ],
    )
    assert np.allclose(
        cc_south.clast__elev.values[:, -1],
        grid_south.at_node["topographic__elevation"][
            cc_south.clast__node.values
        ],
    )

    uplift = np.random.rand(cc_south._grid.number_of_nodes)
    cc_south.run_one_step(
        dt=10.0,
        Si=1.2,
        kappa=kappa,
        uplift=uplift,
        hillslope_river__threshold=1e4,
        lateral_spreading="on",
        disturbance_fqcy=1.0,
        d_star=1.0,
    )
    grid_south.at_node["topographic__elevation"] += (
        10 * uplift
    )

    assert (
        len(cc_south.uplift)
        == cc_south._grid.number_of_nodes
    )
    assert np.allclose(
        cc_south.clast__elev.values[:, -1],
        grid_south.at_node["topographic__elevation"][
            cc_south.clast__node.values
        ],
    )

    cc_south.clast__elev[1, -1] = 0.0
    cc_south.run_one_step(
        dt=10.0,
        Si=1.2,
        kappa=kappa,
        uplift=uplift,
        hillslope_river__threshold=1e4,
        lateral_spreading="on",
        disturbance_fqcy=1.0,
        d_star=1.0,
    )
    grid_south.at_node["topographic__elevation"] += (
        10 * uplift
    )
    assert np.allclose(
        cc_south.clast__elev.values[0, -1],
        grid_south.at_node["topographic__elevation"][
            cc_south.clast__node.values[0]
        ],
    )
    assert cc_south.clast__elev.values[1, -1] < (
        grid_south.at_node["topographic__elevation"][
            cc_south.clast__node.values[1]
        ]
    )

    uplift = [5, 6]
    with pytest.raises(TypeError):
        cc_south.run_one_step(
            dt=10.0,
            Si=1.2,
            kappa=kappa,
            uplift=uplift,
            hillslope_river__threshold=1e4,
            lateral_spreading="on",
            disturbance_fqcy=1.0,
            d_star=1.0,
        )

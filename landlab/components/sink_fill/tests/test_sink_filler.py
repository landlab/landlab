# -*- coding: utf-8 -*-
"""
test_sink_filler:

Created on Tues Oct 20, 2015

@author: dejh
"""
import numpy as np  # for use of np.round
import pytest
from numpy.testing import assert_array_almost_equal, assert_array_equal

from landlab import BAD_INDEX_VALUE as XX, FieldError, RasterModelGrid
from landlab.components import FlowAccumulator, SinkFiller, SinkFillerBarnes


def test_route_to_multiple_error_raised_init():
    mg = RasterModelGrid((10, 10))
    z = mg.add_zeros("node", "topographic__elevation")
    z += mg.x_of_node + mg.y_of_node
    fa = FlowAccumulator(mg, flow_director="MFD")
    fa.run_one_step()
    with pytest.raises(NotImplementedError):
        SinkFillerBarnes(mg)


def test_route_to_multiple_error_raised_run():
    mg = RasterModelGrid((10, 10))
    z = mg.add_zeros("node", "topographic__elevation")
    z += mg.x_of_node + mg.y_of_node
    sfb = SinkFillerBarnes(mg)
    fa = FlowAccumulator(mg, flow_director="MFD")
    fa.run_one_step()
    with pytest.raises(NotImplementedError):
        sfb.run_one_step()


def test_route_to_multiple_error_raised():
    mg = RasterModelGrid((10, 10))
    z = mg.add_zeros("node", "topographic__elevation")
    z += mg.x_of_node + mg.y_of_node
    fa = FlowAccumulator(mg, flow_director="MFD")
    fa.run_one_step()

    with pytest.raises(NotImplementedError):
        SinkFiller(mg)


def test_check_fields(sink_grid1):
    """
    Check to make sure the right fields have been created.
    """
    SinkFillerBarnes(sink_grid1)
    assert_array_equal(
        np.zeros(sink_grid1.number_of_nodes), sink_grid1.at_node["sediment_fill__depth"]
    )
    with pytest.raises(FieldError):
        sink_grid1.at_node["drainage_area"]


def test_get_lake_ext_margin(sink_grid1):
    hf = SinkFiller(sink_grid1)

    lake = np.array([16, 17, 23, 24, 25, 30, 31, 32])
    ext_margin_returned = hf._get_lake_ext_margin(lake)
    ext_margin = np.array(
        [8, 9, 10, 11, 15, 18, 19, 22, 26, 29, 33, 36, 37, 38, 39, 40]
    )
    assert_array_equal(ext_margin_returned, ext_margin)


def test_get_lake_int_margin(sink_grid1):
    hf = SinkFiller(sink_grid1)

    lake = np.array([16, 17, 18, 23, 24, 25, 26, 30, 31, 32])
    ext_margin = np.array(
        [8, 9, 10, 11, 12, 15, 19, 20, 22, 27, 29, 33, 34, 36, 37, 38, 39, 40]
    )
    int_margin_returned = hf._get_lake_int_margin(lake, ext_margin)
    int_margin = np.array([16, 17, 18, 23, 25, 26, 30, 31, 32])
    assert_array_equal(int_margin_returned, int_margin)


def test_drainage_directions_change(sink_grid1):
    hf = SinkFiller(sink_grid1)

    lake = np.array([22, 23])
    old_elevs = np.ones(49, dtype=float)
    old_elevs[lake] = 0.0
    new_elevs = old_elevs.copy()
    new_elevs[40] = 2.0
    cond = hf.drainage_directions_change(lake, old_elevs, new_elevs)
    assert not cond
    new_elevs[23] = 0.5
    cond = hf.drainage_directions_change(lake, old_elevs, new_elevs)
    assert not cond
    new_elevs[23] = 1.0
    cond = hf.drainage_directions_change(lake, old_elevs, new_elevs)
    assert not cond
    new_elevs[23] = 1.2
    cond = hf.drainage_directions_change(lake, old_elevs, new_elevs)
    assert cond


def test_add_slopes(sink_grid1):
    z = sink_grid1.at_node["topographic__elevation"]
    hf = SinkFiller(sink_grid1)

    new_z = z.copy()
    outlet_elev = z[sink_grid1.outlet]
    hf._elev[sink_grid1.lake] = outlet_elev
    rt2 = np.sqrt(2.0)
    slope_to_add = 0.1
    lake_map = np.empty_like(z)
    lake_map.fill(XX)
    lake_map[sink_grid1.lake] = sink_grid1.lake_code
    hf._lf._lake_map = lake_map
    hf.lake_nodes_treated = np.array([], dtype=int)
    dists = sink_grid1.calc_distances_of_nodes_to_point(
        (sink_grid1.node_x[sink_grid1.outlet], sink_grid1.node_y[sink_grid1.outlet])
    )
    new_z[sink_grid1.lake] = outlet_elev
    new_z[sink_grid1.lake] += dists[sink_grid1.lake] * slope_to_add
    # test the ones we can do easily analytically separately
    straight_north = np.array([23, 16])
    off_angle = 24
    elevs_out, lake_out = hf._add_slopes(
        slope_to_add, sink_grid1.outlet, sink_grid1.lake_code
    )
    assert_array_equal(
        slope_to_add * (np.arange(2.0) + 1.0) + outlet_elev, elevs_out[straight_north]
    )
    assert slope_to_add * rt2 + outlet_elev == pytest.approx(elevs_out[off_angle])
    assert_array_equal(new_z, elevs_out)
    assert_array_equal(sink_grid1.lake, lake_out)


def test_filler_flat(sink_grid2):
    """
    Very simple, though possibly degerate, case, filling a 3x3 hole up to
    the flat surface surrounding it.
    """
    hf = SinkFiller(sink_grid2)
    hf.fill_pits()
    assert_array_equal(hf._elev[sink_grid2.lake], np.ones(9, dtype=float))
    assert_array_equal(
        sink_grid2.at_node["topographic__elevation"][sink_grid2.lake],
        np.ones(9, dtype=float),
    )


def test_filler_inclined(sink_grid3):
    """
    Tests a flat fill into an inclined surface, with two holes.
    """
    hf = SinkFiller(sink_grid3)
    hf.fill_pits()
    assert_array_equal(
        sink_grid3.at_node["topographic__elevation"][sink_grid3.lake1],
        np.ones(9, dtype=float) * 4.0,
    )
    assert_array_equal(
        sink_grid3.at_node["topographic__elevation"][sink_grid3.lake2],
        np.ones(4, dtype=float) * 7.0,
    )


def test_filler_inclined2(sink_grid3):
    """
    Tests an inclined fill into an inclined surface, with two holes.
    """
    fr = FlowAccumulator(sink_grid3, flow_director="D8")
    hf = SinkFiller(sink_grid3, apply_slope=True)

    hf.fill_pits()
    hole1 = np.array(
        [
            4.00009091,
            4.00018182,
            4.00027273,
            4.00036364,
            4.00045455,
            4.00054545,
            4.00063636,
            4.00072727,
            4.00081818,
        ]
    )
    hole2 = np.array([7.16666667, 7.33333333, 7.5, 7.66666667])
    assert_array_almost_equal(
        sink_grid3.at_node["topographic__elevation"][sink_grid3.lake1], hole1
    )
    assert_array_almost_equal(
        sink_grid3.at_node["topographic__elevation"][sink_grid3.lake2], hole2
    )
    fr.run_one_step()
    assert sink_grid3.at_node["flow__sink_flag"][sink_grid3.core_nodes].sum() == 0


def test_stupid_shaped_hole(sink_grid4):
    """Tests inclined fill into a surface with a deliberately awkward shape."""
    fr = FlowAccumulator(sink_grid4, flow_director="D8")
    hf = SinkFiller(sink_grid4, apply_slope=True)
    hf.fill_pits()
    hole1 = np.array(
        [
            4.00007692,
            4.00015385,
            4.00023077,
            4.00030769,
            4.00038462,
            4.00046154,
            4.00053846,
            4.00061538,
            4.00069231,
            4.00076923,
            4.00084615,
        ]
    )
    hole2 = np.array([7.4, 7.2, 7.6])

    assert_array_almost_equal(
        sink_grid4.at_node["topographic__elevation"][sink_grid4.lake1], hole1
    )
    assert_array_almost_equal(
        sink_grid4.at_node["topographic__elevation"][sink_grid4.lake2], hole2
    )
    fr.run_one_step()
    assert sink_grid4.at_node["flow__sink_flag"][sink_grid4.core_nodes].sum() == 0


def test_D4_routing(sink_grid5):
    """
    Tests inclined fill into a surface with a deliberately awkward shape.
    This is testing D4 routing.
    """
    fr = FlowAccumulator(sink_grid5, flow_director="D4")
    hf = SinkFiller(sink_grid5, routing="D4", apply_slope=True)
    hf.fill_pits()
    hole1 = np.array(
        [
            4.00016667,
            4.00033333,
            4.0005,
            4.00008333,
            4.00025,
            4.00041667,
            4.000833,
            4.00066667,
            4.00058333,
            4.00075,
            4.334,
        ]
    )
    hole2 = np.array([7.6, 7.2, 7.4])

    assert_array_almost_equal(
        sink_grid5.at_node["topographic__elevation"][sink_grid5.lake1], hole1
    )
    assert_array_almost_equal(
        sink_grid5.at_node["topographic__elevation"][sink_grid5.lake2], hole2
    )
    fr.run_one_step()
    assert sink_grid5.at_node["flow__sink_flag"][sink_grid5.core_nodes].sum() == 0


def test_D4_filling(sink_grid5):
    """
    Tests inclined fill into a surface with a deliberately awkward shape.
    This is testing D4 without inclining the surface.
    """
    hf = SinkFiller(sink_grid5, routing="D4")
    hf.fill_pits()
    hole1 = 4.0 * np.ones_like(sink_grid5.lake1, dtype=float)
    hole1[-1] += 0.001
    hole2 = 7.0 * np.ones_like(sink_grid5.lake2, dtype=float)

    assert_array_almost_equal(
        sink_grid5.at_node["topographic__elevation"][sink_grid5.lake1], hole1
    )
    assert_array_almost_equal(
        sink_grid5.at_node["topographic__elevation"][sink_grid5.lake2], hole2
    )

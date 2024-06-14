#!/usr/env/python

from collections import deque

import numpy as np
import pytest

from landlab import FieldError
from landlab import HexModelGrid
from landlab import RasterModelGrid
from landlab.components import FlowAccumulator
from landlab.components import FlowDirectorDINF
from landlab.components import FlowDirectorSteepest
from landlab.components import LakeMapperBarnes
from landlab.utils import StablePriorityQueue

"""
These tests test specific aspects of LakeMapperBarnes not picked up in the
various docstrings.
"""


def test_route_to_multiple_error_raised_init():
    mg = RasterModelGrid((10, 10))
    z = mg.add_zeros("topographic__elevation", at="node")
    z += mg.x_of_node + mg.y_of_node
    fa = FlowAccumulator(mg, flow_director="MFD")
    fa.run_one_step()
    with pytest.raises(NotImplementedError):
        LakeMapperBarnes(mg)


def test_bad_init_method1():
    rmg = RasterModelGrid((5, 5), xy_spacing=2.0)
    rmg.add_zeros("topographic__elevation", at="node", dtype=float)
    _ = FlowAccumulator(rmg, flow_director="D8")
    with pytest.raises(ValueError):
        LakeMapperBarnes(rmg, method="Nope")


def test_bad_init_method2():
    rmg = RasterModelGrid((5, 5), xy_spacing=2.0)
    rmg.add_zeros("topographic__elevation", at="node", dtype=float)
    _ = FlowAccumulator(rmg, flow_director="D8")
    with pytest.raises(ValueError):
        LakeMapperBarnes(rmg, method="d8")


def test_bad_init_gridmethod():
    hmg = HexModelGrid((30, 29), spacing=3.0)
    hmg.add_zeros("topographic__elevation", at="node", dtype=float)
    _ = FlowAccumulator(hmg, flow_director="Steepest")
    with pytest.raises(ValueError):
        LakeMapperBarnes(hmg, method="D8")


def test_closed_up_grid():
    mg = RasterModelGrid((5, 5))
    for edge in ("left", "right", "top", "bottom"):
        mg.status_at_node[mg.nodes_at_edge(edge)] = mg.BC_NODE_IS_CLOSED
    mg.add_zeros("topographic__elevation", at="node", dtype=float)
    _ = FlowAccumulator(mg, flow_director="D8")
    with pytest.raises(ValueError):
        LakeMapperBarnes(mg)


def test_neighbor_shaping_no_fldir():
    mg = RasterModelGrid((5, 5))
    mg.add_zeros("topographic__elevation", at="node", dtype=float)
    with pytest.raises(FieldError):
        LakeMapperBarnes(mg, method="D8", redirect_flow_steepest_descent=True)


def test_neighbor_shaping_no_creation():
    mg = RasterModelGrid((5, 5))
    mg.add_zeros("topographic__elevation", at="node", dtype=float)
    _ = FlowAccumulator(mg)
    lmb = LakeMapperBarnes(mg, method="D8", redirect_flow_steepest_descent=False)
    with pytest.raises(AttributeError):
        lmb._neighbor_arrays


def test_neighbor_shaping_D8():
    mg = RasterModelGrid((5, 5))
    mg.add_zeros("topographic__elevation", at="node", dtype=float)
    _ = FlowAccumulator(mg)
    lmb = LakeMapperBarnes(mg, method="D8", redirect_flow_steepest_descent=True)
    for arr in (lmb._neighbor_arrays, lmb._link_arrays):
        assert len(arr) == 2
        assert arr[0].shape == (25, 4)
        assert arr[1].shape == (25, 4)
    assert len(lmb._neighbor_lengths) == mg.number_of_d8


def test_neighbor_shaping_D4():
    mg = RasterModelGrid((5, 5))
    mg.add_zeros("topographic__elevation", at="node", dtype=float)
    _ = FlowAccumulator(mg)
    lmb = LakeMapperBarnes(mg, method="Steepest", redirect_flow_steepest_descent=True)
    for arr in (lmb._neighbor_arrays, lmb._link_arrays):
        assert len(arr) == 1
        assert arr[0].shape == (25, 4)
    assert len(lmb._neighbor_lengths) == mg.number_of_links


def test_neighbor_shaping_hex():
    hmg = HexModelGrid((6, 5), spacing=1.0)
    hmg.add_zeros("topographic__elevation", at="node", dtype=float)
    _ = FlowAccumulator(hmg)
    lmb = LakeMapperBarnes(hmg, redirect_flow_steepest_descent=True)
    for arr in (lmb._neighbor_arrays, lmb._link_arrays):
        assert len(arr) == 1
        assert arr[0].shape == (hmg.number_of_nodes, 6)
    assert len(lmb._neighbor_lengths) == hmg.number_of_links


def test_accum_wo_reroute():
    mg = RasterModelGrid((5, 5))
    mg.add_zeros("topographic__elevation", at="node", dtype=float)
    _ = FlowAccumulator(mg)
    with pytest.raises(ValueError):
        LakeMapperBarnes(
            mg,
            method="Steepest",
            redirect_flow_steepest_descent=False,
            reaccumulate_flow=True,
        )


def test_redirect_no_lakes():
    mg = RasterModelGrid((5, 5))
    mg.add_zeros("topographic__elevation", at="node", dtype=float)
    _ = FlowAccumulator(mg)
    with pytest.raises(ValueError):
        LakeMapperBarnes(
            mg, method="D8", track_lakes=False, redirect_flow_steepest_descent=True
        )


def test_route_to_many():
    mg = RasterModelGrid((5, 5))
    mg.add_zeros("topographic__elevation", at="node", dtype=float)
    fd = FlowDirectorDINF(mg, "topographic__elevation")
    _ = FlowAccumulator(mg)
    fd.run_one_step()
    assert mg.at_node["flow__receiver_node"].shape == (mg.number_of_nodes, 2)
    with pytest.raises(NotImplementedError):
        LakeMapperBarnes(mg, method="D8", redirect_flow_steepest_descent=True)


def test_permitted_overfill():
    mg = RasterModelGrid((3, 7))
    for edge in ("top", "right", "bottom"):
        mg.status_at_node[mg.nodes_at_edge(edge)] = mg.BC_NODE_IS_CLOSED
    z = mg.add_zeros("topographic__elevation", at="node", dtype=float)
    z.reshape(mg.shape)[1, 1:-1] = [1.0, 0.2, 0.1, 1.0000000000000004, 1.5]
    _ = FlowAccumulator(mg)
    lmb = LakeMapperBarnes(mg, method="Steepest")
    lmb._closed = mg.zeros("node", dtype=bool)
    lmb._closed[mg.status_at_node == mg.BC_NODE_IS_CLOSED] = True
    open = StablePriorityQueue()
    edges = np.array([7])
    for edgenode in edges:
        open.add_task(edgenode, priority=z[edgenode])
    lmb._closed[edges] = True
    while True:
        try:
            lmb._fill_one_node_to_slant(
                z, mg.adjacent_nodes_at_node, lmb._pit, open, lmb._closed, True
            )
        except KeyError:
            break


def test_no_reroute():
    mg = RasterModelGrid((5, 5), xy_spacing=2.0)
    z = mg.add_zeros("topographic__elevation", at="node", dtype=float)
    z[1] = -1.0
    z[6] = -2.0
    z[19] = -2.0
    z[18] = -1.0
    z[17] = -3.0
    fd = FlowDirectorSteepest(mg)
    fa = FlowAccumulator(mg)
    lmb = LakeMapperBarnes(
        mg,
        method="Steepest",
        fill_flat=True,
        redirect_flow_steepest_descent=True,
        track_lakes=True,
    )
    openq = StablePriorityQueue()

    lake_dict = {1: deque([6]), 18: deque([17])}
    fd.run_one_step()  # fill the director fields
    fa.run_one_step()  # get a drainage_area
    orig_surf = lmb._track_original_surface()
    lmb._redirect_flowdirs(orig_surf, lake_dict, openq)

    assert mg.at_node["flow__receiver_node"][6] == 1
    assert mg.at_node["flow__receiver_node"][17] == 18
    assert mg.at_node["flow__receiver_node"][18] == 19

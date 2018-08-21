#!/usr/env/python

import pytest
import numpy as np

from landlab.components import LakeMapperBarnes
from landlab import RasterModelGrid, HexModelGrid
from landlab import CLOSED_BOUNDARY, FieldError

"""
These tests test specific aspects of LakeMapperBarnes not picked up in the
various docstrings.
"""

def test_bad_init_method1():
    rmg = RasterModelGrid((5, 5), dx=2.)
    rmg.add_zeros('node', 'topographic__elevation', dtype=float)
    with pytest.raises(ValueError):
        lmb = LakeMapperBarnes(rmg, method='Nope')


def test_bad_init_method1():
    rmg = RasterModelGrid((5, 5), dx=2.)
    rmg.add_zeros('node', 'topographic__elevation', dtype=float)
    with pytest.raises(ValueError):
        lmb = LakeMapperBarnes(rmg, method='d8')


def test_bad_init_gridmethod():
    hmg = HexModelGrid(30, 29, dx=3.)
    hmg.add_zeros('node', 'topographic__elevation', dtype=float)
    with pytest.raises(ValueError):
        lmb = LakeMapperBarnes(hmg, method='D8')


def closed_up_grid():
    mg = RasterModelGrid((5, 5), dx=1.)
    for edge in ('left', 'right', 'top', 'bottom'):
        mg.status_at_node[mg.nodes_at_edge(edge)] = CLOSED_BOUNDARY
    mg.add_zeros('node', 'topographic__elevation', dtype=float)
    with pytest.raises(ValueError):
        lmb = LakeMapperBarnes(mg)


def test_neighbor_shaping_no_fldir():
    mg = RasterModelGrid((5, 5), dx=1.)
    mg.add_zeros('node', 'topographic__elevation', dtype=float)
    with pytest.raises(FieldError):
        lmb = LakeMapperBarnes(mg, method='D8',
                               redirect_flow_steepest_descent=True)


def test_neighbor_shaping_no_creation():
    mg = RasterModelGrid((5, 5), dx=1.)
    mg.add_zeros('node', 'topographic__elevation', dtype=float)
    mg.add_zeros('node', 'topographic__steepest_slope', dtype=float)
    mg.add_zeros('node', 'flow__receiver_node', dtype=int)
    mg.add_zeros('node', 'flow__link_to_receiver_node', dtype=int)
    lmb = LakeMapperBarnes(mg, method='D8',
                           redirect_flow_steepest_descent=False)
    with pytest.raises(AttributeError):
        lmb._neighbor_arrays


def test_neighbor_shaping_D8():
    mg = RasterModelGrid((5, 5), dx=1.)
    mg.add_zeros('node', 'topographic__elevation', dtype=float)
    mg.add_zeros('node', 'topographic__steepest_slope', dtype=float)
    mg.add_zeros('node', 'flow__receiver_node', dtype=int)
    mg.add_zeros('node', 'flow__link_to_receiver_node', dtype=int)
    lmb = LakeMapperBarnes(mg, method='D8',
                           redirect_flow_steepest_descent=True)
    for arr in (lmb._neighbor_arrays, lmb._link_arrays):
        assert len(arr) == 2
        assert arr[0].shape == (25, 4)
        assert arr[1].shape == (25, 4)
    assert len(lmb._neighbor_lengths) == mg.number_of_d8


def test_neighbor_shaping_D4():
    mg = RasterModelGrid((5, 5), dx=1.)
    mg.add_zeros('node', 'topographic__elevation', dtype=float)
    mg.add_zeros('node', 'topographic__steepest_slope', dtype=float)
    mg.add_zeros('node', 'flow__receiver_node', dtype=int)
    mg.add_zeros('node', 'flow__link_to_receiver_node', dtype=int)
    lmb = LakeMapperBarnes(mg, method='Steepest',
                           redirect_flow_steepest_descent=True)
    for arr in (lmb._neighbor_arrays, lmb._link_arrays):
        assert len(arr) == 1
        assert arr[0].shape == (25, 4)
    assert len(lmb._neighbor_lengths) == mg.number_of_links


def test_neighbor_shaping_hex():
    hmg = HexModelGrid(6, 5, dx=1.)
    hmg.add_zeros('node', 'topographic__elevation', dtype=float)
    hmg.add_zeros('node', 'topographic__steepest_slope', dtype=float)
    hmg.add_zeros('node', 'flow__receiver_node', dtype=int)
    hmg.add_zeros('node', 'flow__link_to_receiver_node', dtype=int)
    lmb = LakeMapperBarnes(hmg, redirect_flow_steepest_descent=True)
    for arr in (lmb._neighbor_arrays, lmb._link_arrays):
        assert len(arr) == 1
        assert arr[0].shape == (hmg.number_of_nodes, 6)
    assert len(lmb._neighbor_lengths) == hmg.number_of_links


def test_accum_wo_reroute():
    mg = RasterModelGrid((5, 5), dx=1.)
    mg.add_zeros('node', 'topographic__elevation', dtype=float)
    mg.add_zeros('node', 'topographic__steepest_slope', dtype=float)
    mg.add_zeros('node', 'flow__receiver_node', dtype=int)
    mg.add_zeros('node', 'flow__link_to_receiver_node', dtype=int)
    with pytest.raises(ValueError):
        lmb = LakeMapperBarnes(mg, method='Steepest',
                               redirect_flow_steepest_descent=False,
                               reaccumulate_flow=True)

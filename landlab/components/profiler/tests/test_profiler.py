# coding: utf8
# ! /usr/env/python
"""Tests for Profiler.
"""
from matplotlib import pyplot as plt
import numpy as np
import pytest

from landlab import RasterModelGrid
from landlab.components import Profiler


def test_single_segment_profile():
    mg = RasterModelGrid((3, 5))
    core_nodes = mg.core_nodes

    endpoints = [core_nodes.min(), core_nodes.max()]
    profiler = Profiler(mg, endpoints)
    profiler.run_one_step()

    np.testing.assert_array_equal(profiler.network_ids[0], core_nodes)


def test_flipped_single_segment_profile():
    mg = RasterModelGrid((3, 5))
    core_nodes = mg.core_nodes

    endpoints = [core_nodes.max(), core_nodes.min()]
    profiler = Profiler(mg, endpoints)
    profiler.run_one_step()

    np.testing.assert_array_equal(profiler.network_ids[0], np.flip(core_nodes))


def test_steep_profile():
    mg = RasterModelGrid((5, 3))
    core_nodes = mg.core_nodes

    endpoints = [core_nodes.min(), core_nodes.max()]
    profiler = Profiler(mg, endpoints)
    profiler.run_one_step()

    np.testing.assert_array_equal(profiler.network_ids[0], core_nodes)


def test_multi_segment_profile_structure():
    mg = RasterModelGrid((5, 5))
    core_nodes = mg.core_nodes

    endpoints = [core_nodes.min(), core_nodes.min() + 2, core_nodes.max()]
    profiler = Profiler(mg, endpoints)
    profiler.run_one_step()

    np.testing.assert_array_equal(list(profiler.network_structure.keys()),
                                  [0, 1])

    np.testing.assert_array_equal(profiler.network_structure[0]['ids'],
                                  [6, 7, 8])
    np.testing.assert_array_equal(profiler.network_structure[1]['ids'],
                                  [8, 13, 18])

    np.testing.assert_array_equal(profiler.network_structure[0]['distances'],
                                  [0, 1, 2])
    np.testing.assert_array_equal(profiler.network_structure[1]['distances'],
                                  [2, 3, 4])


def test_endpoint_options():
    mg = RasterModelGrid((3, 5))
    core_nodes = mg.core_nodes

    endpoints = [core_nodes.min(), core_nodes.max()]
    profiler_nodes = Profiler(mg, endpoints)
    profiler_nodes.run_one_step()

    ep0 = mg.xy_of_node[core_nodes.min()]
    ep1 = mg.xy_of_node[core_nodes.max()]
    profiler_coords = Profiler(mg, [ep0, ep1])
    profiler_coords.run_one_step()

    np.testing.assert_array_equal(profiler_nodes.network_ids[0],
                                  profiler_coords.network_ids[0])


def test_incorrect_endpoints_type():
    mg = RasterModelGrid((3, 5))
    core_nodes = mg.core_nodes

    endpoints = [core_nodes.min()]
    with pytest.raises(ValueError):
        Profiler(mg, endpoints)


def test_out_of_bounds_point():
    mg = RasterModelGrid((3, 5))
    core_nodes = mg.core_nodes

    endpoints = [core_nodes.min(), 2 * core_nodes.max()]
    with pytest.raises(IndexError):
        Profiler(mg, endpoints)

    ep0 = mg.xy_of_node[core_nodes.min()]
    ep1 = 2 * mg.xy_of_node[core_nodes.max()]
    with pytest.raises(ValueError):
        Profiler(mg, [ep0, ep1])

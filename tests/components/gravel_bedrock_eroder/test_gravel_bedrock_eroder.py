#! /usr/bin/env python
"""
Unit tests for landlab.components.gravel_bedrock_eroder.gravel_bedrock_eroder
"""
import numpy as np
from numpy.testing import assert_almost_equal

from landlab import HexModelGrid, RadialModelGrid, RasterModelGrid
from landlab.components import FlowAccumulator, GravelBedrockEroder


def test_transport_rate():

    grid = HexModelGrid((4, 2), spacing=1000.0)
    grid.status_at_node[grid.perimeter_nodes] = grid.BC_NODE_IS_CLOSED
    grid.status_at_node[0] = grid.BC_NODE_IS_FIXED_VALUE

    elev = grid.add_zeros("topographic__elevation", at="node")
    elev[:] = (0.01 * grid.y_of_node) / np.cos(np.radians(30.0))
    sed = grid.add_zeros("soil__depth", at="node")
    sed[:] = 100.0

    fa = FlowAccumulator(grid)
    fa.run_one_step()
    gbe = GravelBedrockEroder(grid, intermittency_factor=0.02, depth_decay_scale=0.5)
    rock = grid.at_node["bedrock__elevation"]
    qs_out = grid.at_node["bedload_sediment__volume_outflux"]

    gbe.run_one_step(1.0e-6)  # using dt=0 prevents change to sed, rock, or elev
    assert_almost_equal(qs_out[grid.core_nodes], [9.88854526, 3.29618175, 3.29618175])

    elev[:] = (0.01 * grid.y_of_node) / np.cos(np.radians(30.0))
    sed[:] = 0.5
    rock[:] = elev - sed

    gbe.run_one_step(1.0e-6)
    assert_almost_equal(qs_out[grid.core_nodes], [6.25075275, 2.08358425, 2.08358425])

    elev[:] = (0.01 * grid.y_of_node) / np.cos(np.radians(30.0))
    sed[:] = 0.0
    rock[:] = elev

    gbe.run_one_step(1.0e-6)
    assert_almost_equal(qs_out[grid.core_nodes], [0.0, 0.0, 0.0])
    
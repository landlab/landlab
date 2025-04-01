#!/usr/bin/env python3
"""
Created on Thu Apr 23 09:09:49 2020

@author: gtucker
"""

from numpy.testing import assert_array_almost_equal
from numpy.testing import assert_array_equal

from landlab import RasterModelGrid
from landlab.components import FlowAccumulator
from landlab.components import SharedStreamPower


def test_without_depression_handling():
    grid = RasterModelGrid((3, 5), xy_spacing=10.0)
    grid.set_closed_boundaries_at_grid_edges(False, True, False, True)
    z = grid.add_zeros("topographic__elevation", at="node")
    z[grid.x_of_node < 15.0] = 10.0

    fa = FlowAccumulator(grid)
    ed = SharedStreamPower(grid, k_bedrock=0.002, k_transport=0.002)

    fa.run_one_step()
    ed.run_one_step(1.0)

    assert_array_equal(
        ed._q.reshape(grid.shape),
        [
            [0.0, 0.0, 0.0, 0.0, 0.0],
            [0.0, 100.0, 200.0, 100.0, 0.0],
            [0.0, 0.0, 0.0, 0.0, 0.0],
        ],
    )
    assert_array_equal(
        ed._erosion_term.reshape(grid.shape),
        [
            [0.0, 0.0, 0.0, 0.0, 0.0],
            [0.0, 0.02, 0.0, 0.0, 0.0],
            [0.0, 0.0, 0.0, 0.0, 0.0],
        ],
    )
    assert_array_equal(
        ed._depo_rate.reshape(grid.shape),
        [
            [0.0, 0.0, 0.0, 0.0, 0.0],
            [0.0, 0.01, 0.01, 0.0, 0.0],
            [0.0, 0.0, 0.0, 0.0, 0.0],
        ],
    )
    assert_array_equal(
        ed._qs.reshape(grid.shape),
        [
            [0.0, 0.0, 0.0, 0.0, 0.0],
            [0.0, 1.0, 0.0, 0.0, 0.0],
            [0.0, 0.0, 0.0, 0.0, 0.0],
        ],
    )
    assert_array_equal(
        ed.sediment_influx.reshape(grid.shape),
        [
            [0.0, 0.0, 0.0, 0.0, 0.0],
            [0.0, 0.0, 1.0, 0.0, 0.0],
            [0.0, 0.0, 0.0, 0.0, 0.0],
        ],
    )


def test_with_depression_handling():
    grid = RasterModelGrid((3, 5), xy_spacing=10.0)
    grid.set_closed_boundaries_at_grid_edges(False, True, False, True)
    z = grid.add_zeros("topographic__elevation", at="node")
    z[grid.x_of_node < 15.0] = 10.0

    fa = FlowAccumulator(
        grid, routing="D4", depression_finder="DepressionFinderAndRouter"
    )
    ed = SharedStreamPower(grid, k_bedrock=0.002, k_transport=0.002)

    fa.run_one_step()
    ed.run_one_step(1.0)

    assert_array_equal(
        ed._q.reshape(grid.shape),
        [
            [0.0, 0.0, 0.0, 0.0, 0.0],
            [0.0, 100.0, 200.0, 300.0, 300.0],
            [0.0, 0.0, 0.0, 0.0, 0.0],
        ],
    )
    assert_array_equal(
        ed._erosion_term.reshape(grid.shape),
        [
            [0.0, 0.0, 0.0, 0.0, 0.0],
            [0.0, 0.02, 0.0, 0.0, 0.0],
            [0.0, 0.0, 0.0, 0.0, 0.0],
        ],
    )
    assert_array_almost_equal(
        ed._depo_rate.reshape(grid.shape),
        [
            [0.0, 0.0, 0.0, 0.0, 0.0],
            [0.0, 0.01, 0.00333333, 0.00166667, 0.0],
            [0.0, 0.0, 0.0, 0.0, 0.0],
        ],
    )
    assert_array_almost_equal(
        ed._qs.reshape(grid.shape),
        [
            [0.0, 0.0, 0.0, 0.0, 0.0],
            [0.0, 1.0, 0.66666667, 0.5, 0.0],
            [0.0, 0.0, 0.0, 0.0, 0.0],
        ],
    )
    assert_array_almost_equal(
        ed.sediment_influx.reshape(grid.shape),
        [
            [0.0, 0.0, 0.0, 0.0, 0.0],
            [0.0, 0.0, 1.0, 0.66666667, 0.5],
            [0.0, 0.0, 0.0, 0.0, 0.0],
        ],
    )


def test_adaptive_solver_without_depression_handling():
    grid = RasterModelGrid((3, 5), xy_spacing=10.0)
    grid.set_closed_boundaries_at_grid_edges(False, True, False, True)
    z = grid.add_zeros("topographic__elevation", at="node")
    z[grid.x_of_node < 15.0] = 10.0

    fa = FlowAccumulator(grid)
    ed = SharedStreamPower(grid, solver="adaptive", k_bedrock=0.002, k_transport=0.002)

    fa.run_one_step()
    ed.run_one_step(1.0)

    assert_array_equal(
        ed._q.reshape(grid.shape),
        [
            [0.0, 0.0, 0.0, 0.0, 0.0],
            [0.0, 100.0, 200.0, 100.0, 0.0],
            [0.0, 0.0, 0.0, 0.0, 0.0],
        ],
    )
    assert_array_equal(
        ed._erosion_term.reshape(grid.shape),
        [
            [0.0, 0.0, 0.0, 0.0, 0.0],
            [0.0, 0.02, 0.0, 0.0, 0.0],
            [0.0, 0.0, 0.0, 0.0, 0.0],
        ],
    )
    assert_array_equal(
        ed._depo_rate.reshape(grid.shape),
        [
            [0.0, 0.0, 0.0, 0.0, 0.0],
            [0.0, 0.01, 0.01, 0.0, 0.0],
            [0.0, 0.0, 0.0, 0.0, 0.0],
        ],
    )
    assert_array_equal(
        ed._qs.reshape(grid.shape),
        [
            [0.0, 0.0, 0.0, 0.0, 0.0],
            [0.0, 1.0, 0.0, 0.0, 0.0],
            [0.0, 0.0, 0.0, 0.0, 0.0],
        ],
    )
    assert_array_equal(
        ed.sediment_influx.reshape(grid.shape),
        [
            [0.0, 0.0, 0.0, 0.0, 0.0],
            [0.0, 0.0, 1.0, 0.0, 0.0],
            [0.0, 0.0, 0.0, 0.0, 0.0],
        ],
    )


def test_adaptive_solver_with_depression_handling():
    grid = RasterModelGrid((3, 5), xy_spacing=10.0)
    grid.set_closed_boundaries_at_grid_edges(False, True, False, True)
    z = grid.add_zeros("topographic__elevation", at="node")
    z[grid.x_of_node < 15.0] = 10.0

    fa = FlowAccumulator(
        grid, routing="D4", depression_finder="DepressionFinderAndRouter"
    )
    ed = SharedStreamPower(grid, solver="adaptive", k_bedrock=0.002, k_transport=0.002)

    fa.run_one_step()
    ed.run_one_step(1.0)

    assert_array_equal(
        ed._q.reshape(grid.shape),
        [
            [0.0, 0.0, 0.0, 0.0, 0.0],
            [0.0, 100.0, 200.0, 300.0, 300.0],
            [0.0, 0.0, 0.0, 0.0, 0.0],
        ],
    )
    assert_array_equal(
        ed._erosion_term.reshape(grid.shape),
        [
            [0.0, 0.0, 0.0, 0.0, 0.0],
            [0.0, 0.02, 0.0, 0.0, 0.0],
            [0.0, 0.0, 0.0, 0.0, 0.0],
        ],
    )
    assert_array_almost_equal(
        ed._depo_rate.reshape(grid.shape),
        [
            [0.0, 0.0, 0.0, 0.0, 0.0],
            [0.0, 0.01, 0.00333333, 0.00166667, 0.0],
            [0.0, 0.0, 0.0, 0.0, 0.0],
        ],
    )
    assert_array_almost_equal(
        ed._qs.reshape(grid.shape),
        [
            [0.0, 0.0, 0.0, 0.0, 0.0],
            [0.0, 1.0, 0.66666667, 0.5, 0.0],
            [0.0, 0.0, 0.0, 0.0, 0.0],
        ],
    )
    assert_array_almost_equal(
        ed.sediment_influx.reshape(grid.shape),
        [
            [0.0, 0.0, 0.0, 0.0, 0.0],
            [0.0, 0.0, 1.0, 0.66666667, 0.5],
            [0.0, 0.0, 0.0, 0.0, 0.0],
        ],
    )

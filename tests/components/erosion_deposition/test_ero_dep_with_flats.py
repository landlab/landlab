#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Apr 23 09:09:49 2020

@author: gtucker
"""

from numpy.testing import assert_array_almost_equal, assert_array_equal

from landlab import RasterModelGrid
from landlab.components import ErosionDeposition, FlowAccumulator


def test_without_depression_handling():

    grid = RasterModelGrid((3, 5), xy_spacing=10.0)
    grid.set_closed_boundaries_at_grid_edges(False, True, False, True)
    z = grid.add_zeros("node", "topographic__elevation")
    z[grid.x_of_node < 15.0] = 10.0

    fa = FlowAccumulator(grid)
    ed = ErosionDeposition(grid)

    fa.run_one_step()
    ed.run_one_step(1.0)

    assert_array_equal(
        ed._q,
        [
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            100.0,
            200.0,
            100.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
        ],
    )
    assert_array_equal(
        ed._erosion_term,
        [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.02, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],
    )
    assert_array_equal(
        ed._depo_rate,
        [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.01, 0.01, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],
    )
    assert_array_equal(
        ed._qs,
        [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],
    )
    assert_array_equal(
        ed._qs_in,
        [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],
    )


def test_with_depression_handling():

    grid = RasterModelGrid((3, 5), xy_spacing=10.0)
    grid.set_closed_boundaries_at_grid_edges(False, True, False, True)
    z = grid.add_zeros("node", "topographic__elevation")
    z[grid.x_of_node < 15.0] = 10.0

    fa = FlowAccumulator(
        grid, routing="D4", depression_finder="DepressionFinderAndRouter"
    )
    ed = ErosionDeposition(grid)

    fa.run_one_step()
    ed.run_one_step(1.0)

    assert_array_equal(
        ed._q,
        [
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            100.0,
            200.0,
            300.0,
            300.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
        ],
    )
    assert_array_equal(
        ed._erosion_term,
        [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.02, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],
    )
    assert_array_almost_equal(
        ed._depo_rate,
        [
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.01,
            0.00333333,
            0.00166667,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
        ],
    )
    assert_array_almost_equal(
        ed._qs,
        [
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            1.0,
            0.66666667,
            0.5,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
        ],
    )
    assert_array_almost_equal(
        ed._qs_in,
        [
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            1.0,
            0.66666667,
            0.5,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
        ],
    )


def test_adaptive_solver_without_depression_handling():

    grid = RasterModelGrid((3, 5), xy_spacing=10.0)
    grid.set_closed_boundaries_at_grid_edges(False, True, False, True)
    z = grid.add_zeros("node", "topographic__elevation")
    z[grid.x_of_node < 15.0] = 10.0

    fa = FlowAccumulator(grid)
    ed = ErosionDeposition(grid, solver="adaptive")

    fa.run_one_step()
    ed.run_one_step(1.0)

    assert_array_equal(
        ed._q,
        [
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            100.0,
            200.0,
            100.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
        ],
    )
    assert_array_equal(
        ed._erosion_term,
        [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.02, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],
    )
    assert_array_equal(
        ed._depo_rate,
        [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.01, 0.01, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],
    )
    assert_array_equal(
        ed._qs,
        [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],
    )
    assert_array_equal(
        ed._qs_in,
        [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],
    )


def test_adaptive_solver_with_depression_handling():

    grid = RasterModelGrid((3, 5), xy_spacing=10.0)
    grid.set_closed_boundaries_at_grid_edges(False, True, False, True)
    z = grid.add_zeros("node", "topographic__elevation")
    z[grid.x_of_node < 15.0] = 10.0

    fa = FlowAccumulator(
        grid, routing="D4", depression_finder="DepressionFinderAndRouter"
    )
    ed = ErosionDeposition(grid, solver="adaptive")

    fa.run_one_step()
    ed.run_one_step(1.0)

    assert_array_equal(
        ed._q,
        [
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            100.0,
            200.0,
            300.0,
            300.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
        ],
    )
    assert_array_equal(
        ed._erosion_term,
        [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.02, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],
    )
    assert_array_almost_equal(
        ed._depo_rate,
        [
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.01,
            0.00333333,
            0.00166667,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
        ],
    )
    assert_array_almost_equal(
        ed._qs,
        [
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            1.0,
            0.66666667,
            0.5,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
        ],
    )
    assert_array_almost_equal(
        ed._qs_in,
        [
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            1.0,
            0.66666667,
            0.5,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
        ],
    )

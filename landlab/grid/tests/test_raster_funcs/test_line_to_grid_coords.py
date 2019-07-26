#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jul 26 13:59:55 2019

@author: gtucker
"""

from landlab.grid.raster_funcs import line_to_grid_coords
from numpy.testing import assert_array_equal


def test_line_to_grid_coords():
    """Compass directions refer to orientation tested."""
    xy = line_to_grid_coords(2, 0, 5, 0)  # E
    assert_array_equal(xy, [[2, 0],
                            [3, 0],
                            [4, 0],
                            [5, 0]])
    xy = line_to_grid_coords(0, 0, 4, 1)  # ENE
    assert_array_equal(xy, [[0, 0],
                            [1, 0],
                            [2, 0],
                            [3, 1],
                            [4, 1]])
    xy = line_to_grid_coords(0, 2, 1, 5)  # NNE
    assert_array_equal(xy, [[0, 2],
                            [0, 3],
                            [1, 4],
                            [1, 5]])
    xy = line_to_grid_coords(-2, 0, -2, 8)  # N
    assert_array_equal(xy, [[-2, 0],
                            [-2, 1],
                            [-2, 2],
                            [-2, 3],
                            [-2, 4],
                            [-2, 5],
                            [-2, 6],
                            [-2, 7],
                            [-2, 8]])
    xy = line_to_grid_coords(1, 1, 0, 5)  # NNW
    assert_array_equal(xy, [[1, 1],
                            [1, 2],
                            [0, 3],
                            [0, 4],
                            [0, 5]])
    xy = line_to_grid_coords(0, -1, -5, 0)  # WNW
    assert_array_equal(xy, [[0, -1],
                            [-1, -1],
                            [-2, -1],
                            [-3, 0],
                            [-4, 0],
                            [-5, 0]])
    xy = line_to_grid_coords(0, -1, -6, -1)  # W
    assert_array_equal(xy, [[0, -1],
                            [-1, -1],
                            [-2, -1],
                            [-3, -1],
                            [-4, -1],
                            [-5, -1],
                            [-6, -1]])
    xy = line_to_grid_coords(-7, -1, -13, -4)  # WSW
    assert_array_equal(xy, [[-7, -1],
                            [-8, -2],
                            [-9, -2],
                            [-10, -2],
                            [-11, -3],
                            [-12, -4],
                            [-13, -4]])
    xy = line_to_grid_coords(9, -4, 6, -11)  # SSW
    assert_array_equal(xy, [[9, -4],
                            [9, -5],
                            [8, -6],
                            [8, -7],
                            [7, -8],
                            [7, -9],
                            [6, -10],
                            [6, -11]])
    xy = line_to_grid_coords(-2, 1, -2, -6)  # S
    assert_array_equal(xy, [[-2, 1],
                            [-2, 0],
                            [-2, -1],
                            [-2, -2],
                            [-2, -3],
                            [-2, -4],
                            [-2, -5],
                            [-2, -6]])
    xy = line_to_grid_coords(0, 0, 3, -6)  # SSE
    assert_array_equal(xy, [[0, 0],
                            [0, -1],
                            [1, -2],
                            [2, -3],
                            [2, -4],
                            [2, -5],
                            [3, -6]])
    xy = line_to_grid_coords(0, 0, 7, -3)  # ESE
    assert_array_equal(xy, [[0, 0],
                            [1, 0],
                            [2, -1],
                            [3, -1],
                            [4, -2],
                            [5, -2],
                            [6, -3],
                            [7, -3]])

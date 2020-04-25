#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jul 26 13:59:55 2019

@author: gtucker
"""

from numpy.testing import assert_array_equal

from landlab.grid.raster_funcs import line_to_grid_coords


def test_line_to_grid_coords():
    """Compass directions refer to orientation tested."""
    rr, cc = line_to_grid_coords(2, 0, 5, 0)  # E
    assert_array_equal(cc, [2, 3, 4, 5])
    assert_array_equal(rr, [0, 0, 0, 0])

    rr, cc = line_to_grid_coords(0, 0, 4, 1)  # ENE
    assert_array_equal(cc, [0, 1, 2, 3, 4])
    assert_array_equal(rr, [0, 0, 0, 1, 1])

    rr, cc = line_to_grid_coords(0, 2, 1, 5)  # NNE
    assert_array_equal(cc, [0, 0, 1, 1])
    assert_array_equal(rr, [2, 3, 4, 5])

    rr, cc = line_to_grid_coords(-2, 0, -2, 8)  # N
    assert_array_equal(cc, [-2, -2, -2, -2, -2, -2, -2, -2, -2])
    assert_array_equal(rr, [0, 1, 2, 3, 4, 5, 6, 7, 8])

    rr, cc = line_to_grid_coords(1, 1, 0, 5)  # NNW
    assert_array_equal(cc, [1, 1, 0, 0, 0])
    assert_array_equal(rr, [1, 2, 3, 4, 5])

    rr, cc = line_to_grid_coords(0, -1, -5, 0)  # WNW
    assert_array_equal(cc, [0, -1, -2, -3, -4, -5])
    assert_array_equal(rr, [-1, -1, -1, 0, 0, 0])

    rr, cc = line_to_grid_coords(0, -1, -6, -1)  # W
    assert_array_equal(cc, [0, -1, -2, -3, -4, -5, -6])
    assert_array_equal(rr, [-1, -1, -1, -1, -1, -1, -1])

    rr, cc = line_to_grid_coords(-7, -1, -13, -4)  # WSW
    assert_array_equal(cc, [-7, -8, -9, -10, -11, -12, -13])
    assert_array_equal(rr, [-1, -2, -2, -2, -3, -4, -4])

    rr, cc = line_to_grid_coords(9, -4, 6, -11)  # SSW
    assert_array_equal(cc, [9, 9, 8, 8, 7, 7, 6, 6])
    assert_array_equal(rr, [-4, -5, -6, -7, -8, -9, -10, -11])

    rr, cc = line_to_grid_coords(-2, 1, -2, -6)  # S
    assert_array_equal(cc, [-2, -2, -2, -2, -2, -2, -2, -2])
    assert_array_equal(rr, [1, 0, -1, -2, -3, -4, -5, -6])

    rr, cc = line_to_grid_coords(0, 0, 3, -6)  # SSE
    assert_array_equal(cc, [0, 0, 1, 2, 2, 2, 3])
    assert_array_equal(rr, [0, -1, -2, -3, -4, -5, -6])

    rr, cc = line_to_grid_coords(0, 0, 7, -3)  # ESE
    assert_array_equal(cc, [0, 1, 2, 3, 4, 5, 6, 7])
    assert_array_equal(rr, [0, 0, -1, -1, -2, -2, -3, -3])

import numpy as np
from numpy.testing import assert_array_equal
from pytest import approx

from landlab import HexModelGrid


def test_patches_at_link():
    grid = HexModelGrid((3, 2))
    assert_array_equal(
        grid.patches_at_link,
        [
            [0, -1],
            [1, -1],
            [0, 1],
            [0, 2],
            [2, -1],
            [1, 3],
            [2, 4],
            [3, -1],
            [3, 5],
            [4, 5],
            [4, -1],
            [5, -1],
        ],
    )


def test_link_angle():
    grid = HexModelGrid((3, 2))
    assert grid.angle_of_link / np.pi * 3.0 == approx(
        [0.0, 2.0, 1.0, 2.0, 1.0, 0.0, 0.0, 1.0, 2.0, 1.0, 2.0, 0.0], rel=1.0e-5
    )


def test_link_orientation():
    grid = HexModelGrid((3, 3))
    assert_array_equal(
        grid.orientation_of_link,
        [1, 1, 16, 4, 16, 4, 16, 4, 1, 1, 1, 4, 16, 4, 16, 4, 16, 1, 1],
    )
    assert_array_equal(
        grid.orientation_of_link[2],  # re-do for coverage
        16,
    )


def test_parallel_links_at_link():
    grid = HexModelGrid((4, 3))
    target = np.array(
        [
            [-1, -1],
            [-1, -1],
            [-1, -1],
            [-1, 14],
            [-1, 13],
            [-1, 16],
            [-1, 15],
            [-1, -1],
            [-1, 9],
            [8, 10],
            [9, -1],
            [-1, -1],
            [-1, 25],
            [4, 24],
            [3, 27],
            [6, 26],
            [5, 29],
            [-1, 28],
            [-1, -1],
            [-1, 20],
            [19, 21],
            [20, 22],
            [21, -1],
            [-1, -1],
            [13, -1],
            [12, -1],
            [15, -1],
            [14, -1],
            [17, -1],
            [16, -1],
            [-1, -1],
            [-1, -1],
            [-1, -1],
            [-1, -1],
        ]
    )
    assert_array_equal(grid.parallel_links_at_link, target)


def test_parallel_links_at_link_vertical_orientation():
    grid = HexModelGrid((3, 4), orientation="vertical")
    target = np.array(
        [
            [-1, -1],
            [-1, -1],
            [-1, 12],
            [-1, -1],
            [-1, 10],
            [9, -1],
            [-1, 16],
            [-1, -1],
            [-1, 14],
            [13, 5],
            [4, -1],
            [-1, -1],
            [2, 22],
            [-1, 9],
            [8, 20],
            [19, -1],
            [6, 26],
            [-1, -1],
            [-1, 24],
            [23, 15],
            [14, -1],
            [-1, -1],
            [12, 31],
            [-1, 19],
            [18, 30],
            [29, -1],
            [16, -1],
            [-1, -1],
            [-1, -1],
            [-1, 25],
            [24, -1],
            [22, -1],
            [-1, -1],
            [-1, -1],
        ]
    )
    assert_array_equal(grid.parallel_links_at_link, target)

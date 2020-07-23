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

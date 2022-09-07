import numpy as np
from numpy.testing import assert_array_equal
from pytest import approx

from landlab import FramedVoronoiGrid


def test_rect_patches_at_link():
    grid = FramedVoronoiGrid((3, 2), node_layout="rect", random_seed=False)
    assert_array_equal(
        grid.patches_at_link,
        [[0, -1], [0, -1], [0, 1], [1, -1], [1, 2], [2, -1], [2, 3], [3, -1], [3, -1]],
    )


def test_rect_link_angle():
    grid = FramedVoronoiGrid((3, 3), node_layout="rect", random_seed=False)
    assert grid.angle_of_link[0:2] / np.pi == approx(
        np.array([0.0, 0.0, 0.5]),
        rel=1.0e-5,
    )

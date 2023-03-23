import numpy as np
from numpy.testing import assert_array_equal
from pytest import approx, raises

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
    grid = HexModelGrid((3, 2))
    assert_array_equal(
        np.where(grid.link_has_orientation("e"))[0],
        [0, 5, 6, 11]
    )
    assert_array_equal(
        np.where(grid.link_has_orientation("nne"))[0],
        [2, 4, 7, 9]
    )
    assert_array_equal(
        np.where(grid.link_has_orientation("nnw"))[0],
        [1, 3, 8, 10]
    )
    with raises(ValueError):
        grid.link_has_orientation("ene")  # not a valid orientation

    orient = np.zeros(grid.number_of_links, dtype=bool)
    grid = HexModelGrid((2, 3), orientation='vertical')
    assert_array_equal(
        np.where(grid.link_has_orientation("ene", out=orient))[0],
        [1, 3, 8, 10]
    )
    assert_array_equal(
        np.where(grid.link_has_orientation("n"))[0],
        [2, 5, 6, 9]
    )
    assert_array_equal(
        np.where(grid.link_has_orientation("ese"))[0],
        [0, 4, 7, 11]
    )

    wrong_size_out = np.zeros(grid.number_of_nodes, dtype=bool)
    with raises(ValueError):
        grid.link_has_orientation("ene", out=wrong_size_out)
    wrong_type_out = np.zeros(grid.number_of_links)
    with raises(ValueError):
        grid.link_has_orientation("ene", out=wrong_type_out)

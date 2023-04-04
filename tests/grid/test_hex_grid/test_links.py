import numpy as np
import pytest
from numpy.testing import assert_array_equal
from pytest import approx

from landlab import HexModelGrid
from landlab.grid.linkorientation import LinkOrientation


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


ORIENTATION_ANGLE = {
    LinkOrientation.E: 0.0,
    LinkOrientation.ENE: 30.0,
    LinkOrientation.NNE: 60.0,
    LinkOrientation.N: 90.0,
    LinkOrientation.NNW: 120.0,
    LinkOrientation.ESE: 330.0,
}


@pytest.mark.parametrize("orientation", ("horizontal", "vertical"))
@pytest.mark.parametrize("node_layout", ("hex", "rect"))
@pytest.mark.parametrize("link_orientation", list(ORIENTATION_ANGLE))
def test_link_orientation(orientation, node_layout, link_orientation):
    grid = HexModelGrid((5, 5), orientation=orientation, node_layout=node_layout)

    links = grid.orientation_of_link == link_orientation
    assert np.all(
        np.isclose(
            grid.angle_of_link[links],
            ORIENTATION_ANGLE[link_orientation] * np.pi / 180.0,
        )
    )


def test_link_orientation_is_cached():
    """Check successive calls returns cached array."""
    grid = HexModelGrid((3, 3))
    assert grid.orientation_of_link is grid.orientation_of_link


def test_link_orientation_default_hex():
    grid = HexModelGrid((3, 3))
    assert_array_equal(
        grid.orientation_of_link,
        [
            LinkOrientation.E,
            LinkOrientation.E,
            LinkOrientation.NNW,
            LinkOrientation.NNE,
            LinkOrientation.NNW,
            LinkOrientation.NNE,
            LinkOrientation.NNW,
            LinkOrientation.NNE,
            LinkOrientation.E,
            LinkOrientation.E,
            LinkOrientation.E,
            LinkOrientation.NNE,
            LinkOrientation.NNW,
            LinkOrientation.NNE,
            LinkOrientation.NNW,
            LinkOrientation.NNE,
            LinkOrientation.NNW,
            LinkOrientation.E,
            LinkOrientation.E,
        ],
    )
    assert_array_equal(
        grid.orientation_of_link[2],  # re-do for coverage
        LinkOrientation.NNW,
    )

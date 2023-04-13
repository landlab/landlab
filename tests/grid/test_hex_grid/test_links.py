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


@pytest.mark.parametrize("orientation", ("horizontal", "vertical"))
@pytest.mark.parametrize("node_layout", ("hex", "rect"))
def test_link_orientation_vertical(orientation, node_layout):
    triangle_right = LinkOrientation.N | LinkOrientation.ENE | LinkOrientation.ESE
    triangle_up = LinkOrientation.E | LinkOrientation.NNE | LinkOrientation.NNW

    if orientation == "vertical":
        valid_orientations, invalid_orientations = triangle_right, triangle_up
    else:
        valid_orientations, invalid_orientations = triangle_up, triangle_right

    grid = HexModelGrid((5, 5), orientation=orientation, node_layout=node_layout)

    valid_links = grid.orientation_of_link & valid_orientations
    assert np.all(grid.orientation_of_link[valid_links.astype(bool)])

    invalid_links = grid.orientation_of_link & invalid_orientations
    assert not np.any(grid.orientation_of_link[invalid_links.astype(bool)])


def test_link_orientation_is_cached():
    """Check successive calls returns the same cached array."""
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


@pytest.mark.parametrize("orientation", ("horizontal", "vertical"))
@pytest.mark.parametrize("node_layout", ("hex", "rect"))
def test_parallel_links_at_link_same_orientation(orientation, node_layout):
    """Check that parallel links have the same orientation."""
    grid = HexModelGrid((5, 5), orientation=orientation, node_layout=node_layout)

    orientation_of_parallels = grid.orientation_of_link[grid.parallel_links_at_link]
    orientation_of_parallels[grid.parallel_links_at_link == -1] = 0
    out = np.full_like(orientation_of_parallels, True)

    np.equal(
        orientation_of_parallels.T,
        grid.orientation_of_link,
        where=orientation_of_parallels.T > 0,
        out=out.T,
    )

    assert np.all(out)


def test_parallel_links_at_link_horizontal_orientation():
    grid = HexModelGrid((4, 3))

    assert_array_equal(
        grid.parallel_links_at_link[grid.orientation_of_link == LinkOrientation.E],
        [
            [-1, 1],
            [0, -1],
            [-1, 9],
            [8, 10],
            [9, -1],
            [-1, 20],
            [19, 21],
            [20, 22],
            [21, -1],
            [-1, 32],
            [31, 33],
            [32, -1],
        ],
    )
    assert_array_equal(
        grid.parallel_links_at_link[grid.orientation_of_link == LinkOrientation.NNE],
        [
            [-1, 14],
            [-1, 16],
            [-1, 18],
            [-1, 25],
            [3, 27],
            [5, 29],
            [7, -1],
            [-1, -1],
            [12, -1],
            [14, -1],
            [16, -1],
        ],
    )
    assert_array_equal(
        grid.parallel_links_at_link[grid.orientation_of_link == LinkOrientation.NNW],
        [
            [-1, 11],
            [-1, 13],
            [-1, 15],
            [2, -1],
            [4, 24],
            [6, 26],
            [-1, 28],
            [13, -1],
            [15, -1],
            [17, -1],
            [-1, -1],
        ],
    )


def test_parallel_links_at_link_vertical_orientation():
    grid = HexModelGrid((3, 4), orientation="vertical")

    assert_array_equal(
        grid.parallel_links_at_link[grid.orientation_of_link == LinkOrientation.ENE],
        [
            [-1, -1],
            [-1, 10],
            [-1, 14],
            [4, -1],
            [8, 20],
            [-1, 24],
            [14, -1],
            [18, 30],
            [-1, 32],
            [24, -1],
            [28, -1],
        ],
    )

    assert_array_equal(
        grid.parallel_links_at_link[grid.orientation_of_link == LinkOrientation.N],
        [
            [-1, 12],
            [-1, 16],
            [-1, 17],
            [-1, 21],
            [2, 22],
            [6, 26],
            [7, 27],
            [11, -1],
            [12, 31],
            [16, -1],
            [17, -1],
            [22, -1],
        ],
    )

    assert_array_equal(
        grid.parallel_links_at_link[grid.orientation_of_link == LinkOrientation.ESE],
        [
            [3, -1],
            [-1, 0],
            [9, -1],
            [13, 5],
            [-1, 9],
            [19, -1],
            [23, 15],
            [-1, 19],
            [29, -1],
            [-1, 25],
            [-1, -1],
        ],
    )

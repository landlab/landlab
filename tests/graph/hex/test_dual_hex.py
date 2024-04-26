"""Test HexGraph and DualHexGraph."""

import numpy as np
from numpy.testing import assert_array_almost_equal
from numpy.testing import assert_array_equal
from pytest import approx

from landlab.graph import DualHexGraph

ROOT_3_OVER_2 = np.sqrt(3.0) * 0.5


def test_create():
    """Test creating a dual hex graph with rectangular layout."""
    graph = DualHexGraph((4, 3), node_layout="rect")

    assert graph.number_of_nodes == 12
    assert graph.number_of_links == 23
    assert graph.number_of_patches == 12

    assert graph.number_of_corners == 12
    assert graph.number_of_faces == 13
    assert graph.number_of_cells == 2


def test_create_hex():
    """Test creating a dual hex graph with hex layout."""
    graph = DualHexGraph((4, 3), node_layout="hex")

    assert graph.number_of_nodes == 16
    assert graph.number_of_links == 34
    assert graph.number_of_patches == 19

    assert graph.number_of_corners == 19
    assert graph.number_of_faces == 23
    assert graph.number_of_cells == 5


def test_origin():
    """Test setting the origin."""
    graph = DualHexGraph((4, 3))

    assert graph.y_of_node.min() == approx(0.0)
    assert graph.x_of_node.min() == approx(0.0)
    # assert graph.x_of_corner[0] == approx(1.5)

    graph = DualHexGraph((4, 3), xy_of_lower_left=(0.25, 0.5))

    assert graph.y_of_node.min() == approx(0.5)
    assert graph.x_of_node.min() == approx(0.25)
    # assert graph.x_of_corner[0] == approx(1.75)


def test_orientation():
    """Test vertical and horizontal orientation."""
    graph = DualHexGraph((3, 4), orientation="vertical", sort=True)
    assert_array_almost_equal(
        graph.y_of_corner, [0.5, 0.5, 0.5, 1, 1, 1, 1.5, 1.5, 1.5, 2, 2, 2]
    )

    graph = DualHexGraph((4, 3), orientation="horizontal", sort=True)
    assert_array_almost_equal(
        graph.x_of_corner, [0.5, 1.5, 1, 2, 1, 2, 0.5, 1.5, 0.5, 1.5, 1, 2]
    )


def test_adjacent_corners_at_corner():
    graph = DualHexGraph((3, 3), node_layout="hex", sort=True)
    assert_array_equal(
        graph.adjacent_corners_at_corner,
        [
            [3, 2, -1],
            [4, 3, -1],
            [5, 0, -1],
            [6, 0, 1],
            [7, 1, -1],
            [8, 2, -1],
            [9, 8, 3],
            [9, 4, -1],
            [5, 6, -1],
            [6, 7, -1],
        ],
    )

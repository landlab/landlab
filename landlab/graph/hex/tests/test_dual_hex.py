"""Test HexGraph and DualHexGraph."""
from nose.tools import (assert_true, assert_false, assert_equal,
                        assert_almost_equal)
from numpy.testing import assert_array_equal, assert_array_almost_equal
import numpy as np

from landlab.graph import DualHexGraph


ROOT_3_OVER_2 = np.sqrt(3.) * .5


def test_create():
    """Test creating a dual hex graph with rectangular layout."""
    graph = DualHexGraph((4, 3), node_layout='rect')

    assert_equal(graph.number_of_nodes, 12)
    assert_equal(graph.number_of_links, 23)
    assert_equal(graph.number_of_patches, 12)

    assert_equal(graph.number_of_corners, 10)
    assert_equal(graph.number_of_faces, 11)
    assert_equal(graph.number_of_cells, 2)


def test_create_hex():
    """Test creating a dual hex graph with hex layout."""
    graph = DualHexGraph((4, 3), node_layout='hex')

    assert_equal(graph.number_of_nodes, 16)
    assert_equal(graph.number_of_links, 34)
    assert_equal(graph.number_of_patches, 19)

    assert_equal(graph.number_of_corners, 19)
    assert_equal(graph.number_of_faces, 23)
    assert_equal(graph.number_of_cells, 5)


def test_create_rect1():
    """Test creating a dual hex graph."""
    graph = DualHexGraph((4, 3), node_layout='rect1')

    assert_equal(graph.number_of_nodes, 14)
    assert_equal(graph.number_of_links, 28)
    assert_equal(graph.number_of_patches, 15)

    assert_equal(graph.number_of_corners, 13)
    assert_equal(graph.number_of_faces, 15)
    assert_equal(graph.number_of_cells, 3)


def test_spacing():
    """Test spacing of nodes."""
    graph = DualHexGraph((20, 31), node_layout='rect1')
    assert_array_almost_equal(graph.length_of_link, 1.)

    graph = DualHexGraph((31, 20), spacing=2, node_layout='rect1')
    assert_array_almost_equal(graph.length_of_link, 2.)


def test_origin():
    """Test setting the origin."""
    graph = DualHexGraph((4, 3))

    assert_almost_equal(graph.y_of_node[0], 0.)
    assert_almost_equal(graph.x_of_node[0], 0.)
    assert_almost_equal(graph.x_of_corner[0], 1.5)

    graph = DualHexGraph((4, 3), origin=(.5, .25))

    assert_almost_equal(graph.y_of_node[0], .5)
    assert_almost_equal(graph.x_of_node[0], .25)
    assert_almost_equal(graph.x_of_corner[0], 1.75)


def test_orientation():
    """Test vertical and horizontal orientation."""
    graph = DualHexGraph((3, 4), orientation='vertical')
    assert_array_almost_equal(graph.y_of_corner,
                              [.5, .5, 1, 1, 1, 1.5, 1.5, 1.5, 2, 2])

    graph = DualHexGraph((4, 3), orientation='horizontal')
    assert_array_almost_equal(graph.x_of_corner,
                              [1.5, 1, 2, 1, 2, .5, 1.5, .5, 1.5, 1])


def test_adjacent_corners_at_corner():
    graph = DualHexGraph((3, 3), node_layout='hex')
    assert_array_equal(graph.adjacent_corners_at_corner,
                       [[ 3,  2, -1],
                        [ 4 , 3, -1],
                        [ 5,  0, -1],
                        [ 6,  0,  1],
                        [ 7,  1, -1],
                        [ 8,  2, -1],
                        [ 9,  8,  3],
                        [ 9,  4, -1],
                        [ 5,  6, -1],
                        [ 6,  7, -1]])

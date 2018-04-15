"""Test HexGraph and DualHexGraph."""
from nose.tools import (assert_true, assert_false, assert_equal,
                        assert_almost_equal)
from numpy.testing import assert_array_equal, assert_array_almost_equal
from scipy.spatial import Voronoi
import numpy as np

from landlab.graph import HexGraph, DualHexGraph


def test_create_rect():
    """Test creating a hex graph with rectangular layout."""
    graph = HexGraph((3, 2), node_layout='rect')

    assert_equal(graph.number_of_nodes, 6)
    assert_equal(graph.number_of_links, 9)
    assert_equal(graph.number_of_patches, 4)


def test_create_hex():
    """Test creating a hex graph with hex layout."""
    graph = HexGraph((3, 2), node_layout='hex')

    assert_equal(graph.number_of_nodes, 7)
    assert_equal(graph.number_of_links, 12)
    assert_equal(graph.number_of_patches, 6)


def test_create_rect1():
    """Test creating a hex graph."""
    graph = HexGraph((3, 2), node_layout='rect1')

    assert_equal(graph.number_of_nodes, 7)
    assert_equal(graph.number_of_links, 12)
    assert_equal(graph.number_of_patches, 6)


def test_spacing():
    """Test spacing of nodes."""
    graph = HexGraph((20, 31))
    assert_array_almost_equal(graph.length_of_link, 1.)

    graph = HexGraph((31, 20), spacing=2)
    assert_array_almost_equal(graph.length_of_link, 2.)


def test_origin():
    """Test setting the origin."""
    graph = HexGraph((4, 3))

    assert_almost_equal(graph.y_of_node[0], 0.)
    assert_almost_equal(graph.x_of_node[0], 0.)

    graph = HexGraph((4, 3), origin=(.5, .25))

    assert_almost_equal(graph.y_of_node[0], .5)
    assert_almost_equal(graph.x_of_node[0], .25)


def test_orientation():
    """Test vertical and horizontal orientation."""
    graph = HexGraph((3, 3), orientation='vertical')
    assert_array_almost_equal(graph.y_of_node,
                              [0., 0., .5, 1., 1., 1.5, 2., 2., 2.5])

    graph = HexGraph((3, 3), orientation='horizontal')
    assert_array_almost_equal(graph.x_of_node,
                              [.0, 1., 2., .5, 1.5, 2.5, 0., 1., 2.])


def test_perimeter_nodes_rect():
    graph = HexGraph((3, 4), node_layout='rect')
    assert_array_equal(graph.perimeter_nodes, [3, 7, 11, 10, 9, 8, 4, 0, 1, 2])


def test_perimeter_nodes_hex():
    graph = HexGraph((4, 2), node_layout='hex')
    assert_array_equal(graph.perimeter_nodes, [1, 4, 8, 11, 10, 9, 5, 2, 0])


def test_adjacent_nodes_at_node():
    graph = HexGraph((3, 3), node_layout='hex')
    assert_array_equal(graph.adjacent_nodes_at_node,
                       [[ 1,  4,  3, -1, -1, -1],
                        [ 2,  5,  4,  0, -1, -1],
                        [ 6,  5,  1, -1, -1, -1],
                        [ 4,  7,  0, -1, -1, -1],
                        [ 5,  8,  7,  3,  0,  1],
                        [ 6,  9,  8,  4,  1,  2],
                        [ 9,  5,  2, -1, -1, -1],
                        [ 8,  3,  4, -1, -1, -1],
                        [ 9,  7,  4,  5, -1, -1],
                        [ 8,  5,  6, -1, -1, -1]])

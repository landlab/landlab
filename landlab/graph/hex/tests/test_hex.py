"""Test HexGraph and DualHexGraph."""
from nose.tools import (assert_true, assert_false, assert_equal,
                        assert_almost_equal, assert_is_instance)
from numpy.testing import assert_array_equal, assert_array_almost_equal
from scipy.spatial import Voronoi
import numpy as np

from landlab.graph import HexGraph, DualHexGraph


def test_create():
    """Test creating a hex graph."""
    graph = HexGraph((3, 2))

    assert_equal(graph.number_of_nodes, 7)
    assert_equal(graph.number_of_links, 12)
    assert_equal(graph.number_of_patches, 6)


def test_origin():
    """Test spacing of nodes."""
    graph = HexGraph((3, 2), origin=(- np.sqrt(3.) / 2., -1.))

    assert_almost_equal(graph.x_of_node[3], 0.)
    assert_almost_equal(graph.y_of_node[3], 0.)


def test_spacing():
    """Test spacing of nodes."""
    graph = HexGraph((20, 31))
    assert_array_almost_equal(graph.length_of_link, 1.)

    graph = HexGraph((31, 20), spacing=2)
    assert_array_almost_equal(graph.length_of_link, 2.)


def test_orientation():
    """Test vertical and horizontal orientation."""
    graph = HexGraph((2, 3), orientation='vertical')
    assert_array_almost_equal(graph.y_of_node, [0., .5, .5, 1., 1.5, 1.5, 2.])

    graph = HexGraph((3, 2), orientation='horizontal')
    assert_array_almost_equal(graph.x_of_node, [.5, 1.5, 0., 1., 2., .5, 1.5])

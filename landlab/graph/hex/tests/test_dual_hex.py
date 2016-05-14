"""Test HexGraph and DualHexGraph."""
from nose.tools import (assert_true, assert_false, assert_equal,
                        assert_almost_equal)
from numpy.testing import assert_array_equal, assert_array_almost_equal
import numpy as np

from landlab.graph import DualHexGraph


ROOT_3_OVER_2 = np.sqrt(3.) * .5


def test_create():
    """Test creating a dual hex graph."""
    graph = DualHexGraph((3, 2))

    assert_equal(graph.number_of_nodes, 7)
    assert_equal(graph.number_of_links, 12)
    assert_equal(graph.number_of_patches, 6)

    assert_equal(graph.number_of_corners, 6)
    assert_equal(graph.number_of_faces, 6)
    assert_equal(graph.number_of_cells, 1)


def test_origin():
    """Test spacing of nodes."""
    graph = DualHexGraph((3, 2), origin=(- ROOT_3_OVER_2, -1.))

    assert_almost_equal(graph.x_of_node[3], 0.)
    assert_almost_equal(graph.y_of_node[3], 0.)

    assert_array_almost_equal(graph.x_of_corner, [0., -.5, .5, -.5, .5, 0])
    assert_array_almost_equal(graph.y_of_corner / ROOT_3_OVER_2 * 3,
                              [-2., -1., -1., 1., 1., 2.])


def test_spacing():
    """Test spacing of nodes."""
    graph = DualHexGraph((20, 31))
    assert_array_almost_equal(graph.length_of_link, 1.)

    graph = DualHexGraph((31, 20), spacing=2)
    assert_array_almost_equal(graph.length_of_link, 2.)


def test_orientation():
    """Test vertical and horizontal orientation."""
    graph = DualHexGraph((2, 3), orientation='vertical')
    assert_array_almost_equal(graph.y_of_corner, [.5, .5, 1., 1., 1.5, 1.5])

    graph = DualHexGraph((3, 2), orientation='horizontal')
    assert_array_almost_equal(graph.x_of_corner, [1., .5, 1.5, .5, 1.5, 1.])

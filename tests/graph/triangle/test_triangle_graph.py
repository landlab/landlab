"""Tests for the TriangleGraph object."""

import numpy as np
import pytest
from numpy.testing import assert_array_equal

from landlab.graph.triangle.graph import DualTriangleGraph
from landlab.graph.triangle.mesh import TriangleMesh

try:
    TriangleMesh.validate_triangle()
except FileNotFoundError:
    pytestmark = pytest.mark.skip(reason="triangle is not installed")


def test_graph_init():
    graph = DualTriangleGraph(([0, 0, 11, 11], [0, 10, 10, 0]), triangle_opts="pqa1Djevz")

    assert len(graph.x_of_node) == graph.number_of_nodes
    assert len(graph.y_of_node) == graph.number_of_nodes

    assert np.all(graph.x_of_node >= 0.0) and np.all(graph.x_of_node <= 10.0)
    assert np.all(graph.y_of_node >= 0.0) and np.all(graph.y_of_node <= 11.0)
    assert np.all(graph.x_of_corner >= 0.0) and np.all(graph.x_of_corner <= 10.0)
    assert np.all(graph.y_of_corner >= 0.0) and np.all(graph.y_of_corner <= 11.0)

    assert graph.number_of_corners == graph.number_of_patches

    assert np.all(graph.area_of_cell > 0.0)
    assert np.all(graph.area_of_patch > 0.0)
    assert np.all(graph.length_of_link > 0.0)
    assert np.all(graph.length_of_face > 0.0)


def test_graph_init_zero_length_edges():
    with pytest.raises(RuntimeError):
        DualTriangleGraph(([0, 0, 10, 10], [0, 10, 10, 0]), triangle_opts="pqa1Djevz")


def test_raise_error_if_no_interior_nodes(datadir):
    """If no cells are generated, raise a ValueError."""
    with pytest.raises(ValueError):
        DualTriangleGraph.from_shapefile(
            datadir / "polygon_concave.geojson", triangle_opts="pqa100Djevz"
        )


def test_generate_graph_from_geojson(datadir):
    """Test the graph constructor from a geojson file."""
    graph = DualTriangleGraph.from_shapefile(
        datadir / "polygon_concave.geojson", triangle_opts="pqa10Djevz"
    )

    assert graph.number_of_nodes == 25
    assert len(graph.x_of_node) == graph.number_of_nodes
    assert len(graph.y_of_node) == graph.number_of_nodes
    assert graph.number_of_links == 51
    assert graph.number_of_patches == 26
    assert graph.number_of_corners == graph.number_of_patches
    assert graph.number_of_faces == 27
    assert graph.number_of_cells == 1

    assert np.all(graph.area_of_cell > 0.0)
    assert np.all(graph.area_of_patch > 0.0)
    assert np.all(graph.length_of_link > 0.0)
    assert np.all(graph.length_of_face > 0.0)

    assert_array_equal(graph.nodes_at_patch[:3], [[8, 18, 20], [6, 13, 1], [22, 4, 8]])

    assert_array_equal(graph.nodes_at_link[:3], [[8, 18], [18, 20], [20, 8]])

    assert_array_equal(graph.nodes_at_face[:3], [[8, 18], [18, 20], [20, 8]])

    assert_array_equal(graph.corners_at_face[:3], [[0, 18], [0, 5], [1, 16]])

    assert_array_equal(graph.corners_at_cell[0], [22, 8, 24, 21, 25, 9])

    # assert_array_equal(graph.n_corners_at_cell[0], [6])

    assert_array_equal(graph.cell_at_node[24], [0])

    assert_array_equal(graph.links_at_patch[:3], [[0, 1, 2], [3, 4, 5], [6, 7, 8]])

    assert_array_equal(graph.node_at_cell[0], [24])

    assert_array_equal(graph.faces_at_cell[0], [16, 18, 19, 26, 25, 15])


def test_invalid_polygon_error():
    """Invalid geometries should throw an error."""
    ys = [0, 3, 3, 0, 0, 2, 2, 1, 1, 0, 0]
    xs = [0, 0, 3, 3, 2, 2, 1, 1, 2, 2, 0]

    with pytest.raises(ValueError):
        DualTriangleGraph([ys, xs], triangle_opts="pqa0.1Devjz")

"""Tests for the TriangleGraph object."""

import numpy as np
import pytest
from numpy.testing import assert_array_equal

from landlab.graph.triangle import TriangleGraph, TriangleMesh

if not TriangleMesh.validate_triangle():
    pytestmark = pytest.mark.skip(reason="triangle is not installed")


ys = [0, 0, 10, 10]
xs = [0, 10, 10, 0]


def test_graph_init():
    graph = TriangleGraph(np.array([ys, xs]), triangle_opts="pqa1Djevz")

    assert graph.number_of_nodes == 89
    assert len(graph.x_of_node) == graph.number_of_nodes
    assert len(graph.y_of_node) == graph.number_of_nodes
    assert graph.number_of_links == 232
    assert graph.number_of_patches == 144
    assert graph.number_of_corners == graph.number_of_patches
    assert graph.number_of_faces == 200
    assert graph.number_of_cells == 57

    assert_array_equal(
        graph.nodes_at_patch[:3], [[41, 43, 15], [23, 28, 31], [4, 51, 50]]
    )

    assert_array_equal(graph.nodes_at_link[:3], [[41, 43], [43, 15], [15, 41]])

    assert_array_equal(graph.nodes_at_face[:3], [[41, 43], [43, 15], [15, 41]])

    assert_array_equal(graph.corners_at_face[:3], [[0, 66], [0, 65], [0, 6]])

    assert_array_equal(
        graph.corners_at_cell[:3],
        [
            [2, 62, 68, 70, 71, 77, 79, -1, -1],
            [22, 23, 27, 28, 63, 65, -1, -1, -1],
            [88, 89, 95, 101, 115, 136, 138, -1, -1],
        ],
    )

    assert_array_equal(graph.n_corners_at_cell[:3], [9, 9, 9])

    assert_array_equal(graph.cell_at_node[:5], [-1, -1, -1, -1, 0])

    assert_array_equal(graph.links_at_patch[:3], [[0, 1, 2], [3, 4, 5], [6, 7, 8]])

    assert_array_equal(graph.node_at_cell[:3], [4, 9, 10])

    assert_array_equal(
        graph.faces_at_cell[:3],
        [
            [6, 8, 109, 110, 116, 117, 119, -1, -1],
            [54, 55, 56, 62, 63, 112, -1, -1, -1],
            [145, 146, 148, 149, 164, 166, 180, -1, -1],
        ],
    )

    assert_array_equal(
        graph.perimeter_nodes,
        [0,  1,  2,  3,  5,  6,  7,  8, 
         13, 14, 18, 23, 24, 25, 30, 33, 
         36, 37, 39, 57, 59, 63, 66, 68, 
         70, 72, 75, 77, 82, 83, 87, 88]
    )


def test_generate_graph_from_geojson(datadir):
    """Test the graph constructor from a geojson file."""
    graph = TriangleGraph.from_shapefile(
        datadir / "example_geojson_for_graph.geojson", triangle_opts="pqDjevz"
    )

    assert graph.number_of_nodes == 697
    assert len(graph.x_of_node) == graph.number_of_nodes
    assert len(graph.y_of_node) == graph.number_of_nodes
    assert graph.number_of_links == 1396
    assert graph.number_of_patches == 709
    assert graph.number_of_corners == graph.number_of_patches
    assert graph.number_of_faces == 731
    assert graph.number_of_cells == 32

    assert_array_equal(
        graph.nodes_at_patch[:3], [[476, 474, 475], [609, 478, 601], [476, 473, 474]]
    )

    assert_array_equal(graph.nodes_at_link[:3], [[476, 474], [474, 475], [475, 476]])

    assert_array_equal(graph.nodes_at_face[:3], [[476, 474], [474, 475], [475, 476]])

    assert_array_equal(graph.corners_at_face[:3], [[0, 2], [1, 7], [1, 336]])

    assert_array_equal(
        graph.corners_at_cell[:3],
        [
            [224, 227, 229, 239, 293, 320, 412, -1, -1, -1],
            [267, 269, 281, 283, 392, -1, -1, -1, -1, -1],
            [393, 433, 535, 536, 538, 539, -1, -1, -1, -1],
        ],
    )

    assert_array_equal(graph.n_corners_at_cell[:3], [10, 10, 10])

    assert_array_equal(graph.cell_at_node[:3], [-1, -1, -1])

    assert_array_equal(graph.links_at_patch[:3], [[0, 1, 2], [3, 4, 5], [0, 6, 7]])

    assert_array_equal(graph.node_at_cell[:3], [622, 629, 646])

    assert_array_equal(
        graph.faces_at_cell[:3],
        [
            [274, 276, 278, 279, 282, 294, 296, -1, -1, -1],
            [327, 328, 329, 342, 345, -1, -1, -1, -1, -1],
            [484, 485, 531, 612, 613, 614, -1, -1, -1, -1],
        ],
    )


def test_invalid_polygon_error():
    """Invalid geometries should throw an error."""
    ys = [0, 3, 3, 0, 0, 2, 2, 1, 1, 0, 0]
    xs = [0, 0, 3, 3, 2, 2, 1, 1, 2, 2, 0]

    with pytest.raises(ValueError):
        TriangleGraph([ys, xs], triangle_opts="pqa0.1Devjz")

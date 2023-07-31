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
        [
            0,
            1,
            2,
            3,
            5,
            6,
            7,
            8,
            13,
            14,
            18,
            23,
            24,
            25,
            30,
            33,
            36,
            37,
            39,
            57,
            59,
            63,
            66,
            68,
            70,
            72,
            75,
            77,
            82,
            83,
            87,
            88,
        ],
    )


def test_raise_error_if_no_interior_nodes(datadir):
    """If no cells are generated, raise a ValueError."""
    with pytest.raises(ValueError):
        TriangleGraph.from_shapefile(
            datadir / "polygon_concave.geojson", triangle_opts="pqa100Djevz"
        )


def test_generate_graph_from_geojson(datadir):
    """Test the graph constructor from a geojson file."""
    graph = TriangleGraph.from_shapefile(
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

    assert_array_equal(graph.nodes_at_patch[:3], [[8, 18, 20], [13, 1, 6], [4, 8, 22]])

    assert_array_equal(graph.nodes_at_link[:3], [[8, 18], [18, 20], [20, 8]])

    assert_array_equal(graph.nodes_at_face[:3], [[8, 18], [18, 20], [20, 8]])

    assert_array_equal(graph.corners_at_face[:3], [[0, 18], [0, 5], [1, 16]])

    assert_array_equal(graph.corners_at_cell[0], [8, 9, 21, 22, 24, 25])

    assert_array_equal(graph.n_corners_at_cell[0], [6])

    assert_array_equal(graph.cell_at_node[24], [0])

    assert_array_equal(graph.links_at_patch[:3], [[0, 1, 2], [3, 4, 5], [6, 7, 8]])

    assert_array_equal(graph.node_at_cell[0], [24])

    assert_array_equal(graph.faces_at_cell[0], [15, 16, 18, 19, 25, 26])


def test_invalid_polygon_error():
    """Invalid geometries should throw an error."""
    ys = [0, 3, 3, 0, 0, 2, 2, 1, 1, 0, 0]
    xs = [0, 0, 3, 3, 2, 2, 1, 1, 2, 2, 0]

    with pytest.raises(ValueError):
        TriangleGraph([ys, xs], triangle_opts="pqa0.1Devjz")

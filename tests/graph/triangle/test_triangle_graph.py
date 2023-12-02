"""Tests for the TriangleGraph object."""

import numpy as np
import pytest
from numpy.testing import assert_array_equal

from landlab.graph.triangle.dual_triangle import DualTriGraph
from landlab.graph.triangle.triangle_mesh import TriangleMesh

if not TriangleMesh.validate_triangle():
    pytestmark = pytest.mark.skip(reason="triangle is not installed")


def test_graph_init():
    graph = DualTriGraph(([0, 0, 10, 10], [0, 10, 10, 0]), triangle_opts="pqa1Djevz")

    assert graph.number_of_nodes == 89
    assert len(graph.x_of_node) == graph.number_of_nodes
    assert len(graph.y_of_node) == graph.number_of_nodes
    assert graph.number_of_links == 232
    assert graph.number_of_patches == 144
    assert graph.number_of_corners == graph.number_of_patches
    assert graph.number_of_faces == 200
    assert graph.number_of_cells == 57

    assert np.all(graph.area_of_cell > 0.0)
    assert np.all(graph.area_of_patch > 0.0)
    assert np.all(graph.length_of_link > 0.0)
    assert np.all(graph.length_of_face > 0.0)

    assert_array_equal(
        graph.nodes_at_patch[:3], [[41, 43, 15], [28, 31, 23], [51, 50, 4]]
    )

    assert_array_equal(graph.nodes_at_link[:3], [[41, 43], [43, 15], [15, 41]])

    assert_array_equal(graph.nodes_at_face[:3], [[41, 43], [43, 15], [15, 41]])

    assert_array_equal(graph.corners_at_face[:3], [[0, 66], [0, 65], [0, 6]])

    assert_array_equal(
        graph.corners_at_cell[:3],
        [
            [70, 79, 2, 77, 71, 68, 62, -1, -1],
            [23, 65, 63, 27, 22, 28, -1, -1, -1],
            [88, 115, 136, 101, 138, 89, 95, -1, -1],
        ],
    )

    # assert_array_equal(graph.n_corners_at_cell[:3], [9, 9, 9])

    assert_array_equal(graph.cell_at_node[:5], [-1, -1, -1, -1, 0])

    assert_array_equal(graph.links_at_patch[:3], [[0, 1, 2], [3, 4, 5], [6, 7, 8]])

    assert_array_equal(graph.node_at_cell[:3], [4, 9, 10])

    assert_array_equal(
        graph.faces_at_cell[:3],
        [
            [117, 109, 110, 116, 119, 8, 6, -1, -1],
            [56, 55, 54, 62, 63, 112, -1, -1, -1],
            [180, 145, 146, 148, 149, 166, 164, -1, -1],
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
        DualTriGraph.from_shapefile(
            datadir / "polygon_concave.geojson", triangle_opts="pqa100Djevz"
        )


def test_generate_graph_from_geojson(datadir):
    """Test the graph constructor from a geojson file."""
    graph = DualTriGraph.from_shapefile(
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
        DualTriGraph([ys, xs], triangle_opts="pqa0.1Devjz")

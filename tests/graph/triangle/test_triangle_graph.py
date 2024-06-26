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


@pytest.fixture
def square_graph():
    return DualTriangleGraph(
        ([0, 0, 11, 11], [0, 10, 10, 0]), triangle_opts="pqa1Djevz"
    )


@pytest.mark.parametrize("point", ("corner", "node"))
def test_all_points_in_box(square_graph, point):
    x, y = getattr(square_graph, f"x_of_{point}"), getattr(
        square_graph, f"y_of_{point}"
    )

    assert np.all(x >= 0.0) and np.all(x <= 10.0)
    assert np.all(y >= -1.0) and np.all(y <= 11.0)


@pytest.mark.parametrize("edge", ("face", "link"))
def test_no_zero_length_edges(square_graph, edge):
    assert np.all(getattr(square_graph, f"length_of_{edge}") >= 0.0)


@pytest.mark.parametrize("polygon", ("cell", "patch"))
def test_no_zero_area_polygons(square_graph, polygon):
    assert np.all(getattr(square_graph, f"area_of_{polygon}") >= 0.0)


def test_graph_init():
    graph = DualTriangleGraph(
        ([0, 0, 11, 11], [0, 10, 10, 0]), triangle_opts="pqa1Djevz"
    )
    assert graph.number_of_corners == graph.number_of_patches


def test_graph_init_zero_length_edges():
    with pytest.raises(RuntimeError):
        DualTriangleGraph(([0, 0, 10, 10], [0, 10, 10, 0]), triangle_opts="pqa1Djevz")


def test_raise_error_if_no_interior_nodes(geojson_concave_polygon):
    """If no cells are generated, raise a ValueError."""
    with pytest.raises(ValueError):
        DualTriangleGraph.from_shapefile(
            geojson_concave_polygon, triangle_opts="pqa100Djevz"
        )


def test_circular_polygon(geojson_circular_polygon):
    graph = DualTriangleGraph.from_shapefile(
        geojson_circular_polygon,
        triangle_opts="pqDjevz",
    )

    assert np.all(graph.area_of_cell > 0.0)
    assert np.all(graph.area_of_patch > 0.0)
    assert np.all(graph.length_of_link > 0.0)
    assert np.all(graph.length_of_face > 0.0)

    assert graph.number_of_nodes == 104
    assert graph.number_of_links == 245
    assert graph.number_of_patches == 142
    assert graph.number_of_corners == graph.number_of_patches
    assert graph.number_of_faces == 181
    assert graph.number_of_cells == 40


def test_generate_graph_from_geojson(geojson_concave_polygon):
    """Test the graph constructor from a geojson file."""
    graph = DualTriangleGraph.from_shapefile(
        geojson_concave_polygon, triangle_opts="pqa10Djevz"
    )

    assert graph.number_of_nodes == 25
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


def test_multiple_interior_rings(geojson_interior_rings):
    graph = DualTriangleGraph.from_shapefile(
        geojson_interior_rings, triangle_opts="pq10a10Djevz"
    )

    assert np.all(graph.x_of_node >= 0.0) and np.all(graph.x_of_node <= 10.0)
    assert np.all(graph.y_of_node >= 0.0) and np.all(graph.y_of_node <= 10.0)

    assert not np.any(
        (graph.x_of_node > 3.0)
        & (graph.x_of_node < 4.0)
        & (graph.y_of_node > 3.0)
        & (graph.y_of_node < 4.0)
    )
    assert not np.any(
        (graph.x_of_node > 6.0)
        & (graph.x_of_node < 7.0)
        & (graph.y_of_node > 6.0)
        & (graph.y_of_node < 7.0)
    )

    assert np.all(graph.length_of_link >= 0.0)
    assert np.all(graph.length_of_face >= 0.0)
    assert np.all(graph.area_of_cell >= 0.0)
    assert np.all(graph.area_of_patch >= 0.0)

    assert graph.number_of_corners == graph.number_of_patches


def test_invalid_polygon_error():
    """Invalid geometries should throw an error."""
    ys = [0, 3, 3, 0, 0, 2, 2, 1, 1, 0, 0]
    xs = [0, 0, 3, 3, 2, 2, 1, 1, 2, 2, 0]

    with pytest.raises(ValueError):
        DualTriangleGraph([ys, xs], triangle_opts="pqa0.1Devjz")

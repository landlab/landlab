"""Tests for the TriangleGraph object."""

import numpy as np
from numpy.testing import assert_array_equal
import pytest
from landlab.graph.triangle import TriangleMesh, TriangleGraph

xy_points = np.array([
    [0.0, 0.0],
    [1.0, 0.0],
    [2.0, 0.0],
    [0.5, 1.0],
    [1.5, 1.0],
    [2.5, 1.0],
    [0.0, 2.0],
    [1.0, 2.0],
    [2.0, 2.0],
    [0.0, 3.0],
    [1.0, 3.0],
    [2.0, 3.0],
    [0.0, 0.0]
])

@pytest.fixture
def mesh_from_points():
    mesh = TriangleMesh.from_points(xy_points, opts="pqa0.1Devz")
    mesh.triangulate()
    return mesh

@pytest.fixture
def graph(mesh_from_points):
    mesh = mesh_from_points
    graph = TriangleGraph(mesh.delaunay, mesh.voronoi)
    return graph

def test_graph_init(graph):
    """Test initialization of the TriangleGraph."""
    assert graph.number_of_nodes == 57
    assert len(graph.x_of_node) == graph.number_of_nodes
    assert len(graph.y_of_node) == graph.number_of_nodes
    assert graph.number_of_links == 108
    assert graph.number_of_patches == 53
    assert graph.number_of_corners == graph.number_of_patches
    assert graph.number_of_faces == 51
    assert graph.number_of_cells == 5

    assert_array_equal(
        graph.nodes_at_patch[:3],
        [[53, 1, 54],
         [48, 18, 54],
         [33, 16, 30]]
    )

    assert_array_equal(
        graph.nodes_at_link[:3],
        [[53, 1],
         [1, 54],
         [54, 53]]
    )

    assert_array_equal(
        graph.nodes_at_face[:3],
        [[53, 1],
         [1, 54],
         [54, 53]]
    )

    assert_array_equal(
        graph.corners_at_face[:3],
        [[0, 47],
         [0, 49],
         [0, 51]]
    )

    # TODO FIX THIS
    # assert_array_equal(
    #     graph.corners_at_cell[:3],
    #     []
    # )

    assert_array_equal(
        graph.n_corners_at_cell[:3],
        [4, 4, 4]
    )

    assert_array_equal(
        graph.cell_at_node[:3],
        [0, -1, -1]
    )

    assert_array_equal(
        graph.links_at_patch[:3],
        [[0, 1, 2],
         [3, 4, 5],
         [6, 7, 8]]
    )

    assert_array_equal(
        graph.node_at_cell[:3],
        [0, 45, 46]
    )

    # TODO FIX THIS
    # assert_array_equal(
    #     graph.faces_at_cell[:3],
    #     []
    # )

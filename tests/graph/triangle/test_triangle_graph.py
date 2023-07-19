"""Tests for the TriangleGraph object."""

import numpy as np
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

def test_graph_init(graph):
    """Test initialization of the TriangleGraph."""
    assert True
"""Tests for the TriangleGraph object."""

import pytest

from landlab.graph.triangle import TriangleMesh, TriangleGraph


@pytest.fixture
def triangle_graph():
    path = "tests/graph/triangle/test_triangle_mesh/example_shapefile.geojson"

    # a10000 gives us 18 nodes in the mesh, rather than 5
    mesh = TriangleMesh.from_shapefile(path, opts="pqa10000Dez")
    graph = TriangleGraph(mesh.delaunay, mesh.voronoi)
    return graph


def test_triangle_graph_init(triangle_graph):
    """This test should always pass."""
    assert True

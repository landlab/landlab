"""Tests for the TriangleMesh object."""

import numpy as np
import pytest
from numpy.testing import assert_array_equal, assert_array_almost_equal

from landlab.graph.triangle import TriangleMesh

xy_points = [
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
]

@pytest.fixture
def mesh_from_points():
    mesh = TriangleMesh(xy_points)
    return mesh

@pytest.fixture
def mesh_from_dims():
    shape = (4, 3)
    spacing = (6.0, 2.0)
    mesh = TriangleMesh.from_dims(shape, spacing)
    return mesh

@pytest.fixture
def mesh_from_shapefile():
    path = 'tests/graph/triangle/test_triangle_mesh/test_shapefile.geojson'
    mesh = TriangleMesh.from_shapefile(path)
    return mesh

def test_init_from_points(mesh_from_points):
    """Initialize TriangleMesh from a list of points."""
    mesh = mesh_from_points
    
    assert_array_equal(mesh._vertices, xy_points)

    assert_array_equal(
        mesh._segments,
        [
            [i, i + 1] if i != 11 else [11, 0] 
            for i in range(len(xy_points))
        ]
    )

def test_init_from_dims(mesh_from_dims):
    """Initialize the TriangleMesh from dimension lengths and spacing."""
    mesh = mesh_from_dims

    assert_array_equal(
        mesh._vertices,
        [
            [0, 0],
            [8, 0],
            [16, 0],
            [24, 0],
            [0, 3],
            [8, 3],
            [16, 3],
            [24, 3],
            [0, 6],
            [8, 6],
            [16, 6],
            [24, 6],
        ]
    )

    assert_array_equal(
        mesh._segments,
        [
            [i, i + 1] if i != 11 else [11, 0] 
            for i in range(12)
        ]
    )

def test_init_from_shapefile(mesh_from_shapefile):
    mesh = mesh_from_shapefile

    assert_array_equal(
        mesh._vertices,
        [
            [ 393, -211],
            [ 792, -208],
            [ 803, -530],
            [ 390, -523],
            [ 393, -211]
        ]
    )

    assert_array_equal(
        mesh._segments,
        [
            [i, i + 1] if i != 3 else [3, 0] 
            for i in range(4)
        ]
    )

for mesh in [mesh_from_points, mesh_from_dims, mesh_from_shapefile]:

    def test_triangulate(mesh):

        for required in [
            'vertices', 
            'vertex_markers', 
            'edges',
            'edge_markers',
            'triangles'
        ]:
            assert required in mesh.delaunay.keys()

        for required in [
            'vertices',
            'edges'
        ]:
            assert required in mesh.voronoi.keys()

"""Tests for the TriangleMesh object."""

import numpy as np
import geopandas as gpd

import pytest
from numpy.testing import assert_array_equal

from landlab.graph.triangle import TriangleMesh

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
])

@pytest.fixture
def mesh_from_points():
    mesh = TriangleMesh.from_points(xy_points, opts="pqDez")
    return mesh

@pytest.fixture
def mesh_from_shapefile():
    path = "tests/graph/triangle/test_triangle_mesh/example_geojson.geojson"
    mesh = TriangleMesh.from_shapefile(path, opts="pqDez")
    return mesh

def test_init_from_points(mesh_from_points):
    """Test initialization from list of points."""
    mesh = mesh_from_points

def test_init_from_geojson(mesh_from_shapefile):
    """Test initialization from a geojson file."""
    mesh = mesh_from_shapefile


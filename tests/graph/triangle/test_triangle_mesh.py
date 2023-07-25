"""Tests for the TriangleMesh object."""

import numpy as np
import pytest

from landlab.graph.triangle import TriangleMesh

if not TriangleMesh.validate_triangle():
    pytestmark = pytest.mark.skip(reason="triangle is not installed")


xy_points = np.array(
    [
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
        [0.0, 0.0],
    ]
)


def test_init_from_points():
    """Test initialization from list of points."""
    mesh = TriangleMesh.from_points(xy_points, opts="pqDevjz")

    assert mesh._vertices.shape == (xy_points.shape[0], 2)
    assert mesh._segments.shape == (xy_points.shape[0] - 1, 2)
    assert mesh._holes is None
    assert mesh._opts == "pqDevjz"


def test_triangulate_from_points():
    """Test triangulation routine."""
    mesh = TriangleMesh.from_points(xy_points, opts="pqDevjz")
    mesh.triangulate()


def test_init_from_geojson(datadir):
    """Test initialization from a geojson file."""
    mesh = TriangleMesh.from_shapefile(
        datadir / "example_geojson.geojson", opts="pqDevjz"
    )

    # Validated against Shapely directly
    assert mesh._vertices.shape == (2402, 2)
    assert len(mesh._poly.interiors) == 42
    assert mesh._segments.shape == (2359, 2)
    assert mesh._holes.shape == (42, 2)
    assert mesh._opts == "pqDevjz"


def test_triangulate_from_geojson(datadir):
    """Test triangulation routine."""
    mesh = TriangleMesh.from_shapefile(
        datadir / "example_geojson.geojson", opts="pqDevjz"
    )
    mesh.triangulate()

"""Tests for the TriangleMesh object."""

import numpy as np
import pytest
import shapely
from numpy.testing import assert_array_equal

from landlab.graph.triangle.mesh import TriangleMesh

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


def test_init_from_geojson(concave_polygon):
    """Test initialization from a geojson file."""
    mesh = TriangleMesh.from_shapefile(concave_polygon, opts="pqDevjz")

    assert_array_equal(
        mesh._vertices,
        [
            [0.0, 0.0],
            [10.0, 0.0],
            [10.0, 10.0],
            [5.0, 15.0],
            [0.0, 10.0],
            [0.0, 0.0],
            [2.0, 2.0],
            [8.0, 2.0],
            [8.0, 8.0],
            [2.0, 8.0],
            [2.0, 2.0],
        ],
    )

    assert_array_equal(mesh._holes, [[5.0, 5.0]])

    assert_array_equal(
        mesh._segments,
        [[6, 7], [7, 8], [8, 9], [9, 6], [0, 1], [1, 2], [2, 3], [3, 4], [4, 0]],
    )

    assert mesh._opts == "pqDevjz"


def test_triangulate_from_geojson(concave_polygon):
    """Test triangulation routine."""
    mesh = TriangleMesh.from_shapefile(concave_polygon, opts="pqDevjz")
    mesh.triangulate()


def test_segment(concave_polygon):
    "Test segmentation routine."
    mesh = TriangleMesh.from_shapefile(concave_polygon, opts="pqDevjz")
    segments = mesh._segment(mesh._poly)

    assert len(mesh._holes) == len(mesh._poly.interiors)
    assert len(segments) == 9

    for hole in mesh._holes:
        point = shapely.Point(hole)

        assert not mesh._poly.contains(point)

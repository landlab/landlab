import numpy as np
import pytest
from numpy.testing import assert_array_equal

from landlab.grid.nodestatus import NodeStatus
from landlab.grid.quad import QuadGridGeometry
from landlab.grid.quad import QuadModelGrid
from landlab.grid.quad import _classify_quad_grid_geometry
from landlab.grid.quad import _make_node_status
from landlab.grid.quad import _normalize_bc_at_edge


def test_classify_raster():
    x = [
        [0, 1, 2],
        [0, 1, 2],
        [0, 1, 2],
    ]
    y = [
        [0, 0, 0],
        [2, 2, 2],
        [4, 4, 4],
    ]
    assert _classify_quad_grid_geometry(x, y) == QuadGridGeometry.RASTER


def test_classify_rectlinear():
    x = [
        [0, 1, 3],
        [0, 1, 3],
        [0, 1, 3],
    ]
    y = [
        [0, 0, 0],
        [2, 2, 2],
        [4, 4, 4],
    ]

    assert (
        _classify_quad_grid_geometry(*np.meshgrid(x, y)) == QuadGridGeometry.RECTILINEAR
    )


def test_classify_general():
    x = [
        [0, 1, 2],
        [1, 2, 3],
        [2, 3, 4],
    ]
    y = [
        [0, 0, 0],
        [1, 1, 1],
        [2, 2, 2],
    ]
    assert _classify_quad_grid_geometry(x, y) == QuadGridGeometry.GENERAL


def test_bad_shape():
    with pytest.raises(ValueError, match="coordinate arrays must be 2D"):
        QuadModelGrid(np.arange(10), np.arange(10))


def test_bad_shape_from_rectilinear():
    with pytest.raises(ValueError, match="coordinate arrays must be 1D"):
        QuadModelGrid.from_rectilinear(*np.meshgrid(np.arange(10), np.arange(10)))


def test_shape_mismatch():
    with pytest.raises(ValueError, match="coordinate arrays must have the same shape"):
        QuadModelGrid(np.ones((4, 5)), np.ones((5, 4)))


def test_from_raster():
    grid = QuadModelGrid.from_raster((4, 5), xy_spacing=(1.0, 1.0))

    assert grid.shape == (4, 5)
    assert grid.geometry is QuadGridGeometry.RASTER


def test_from_general():
    x = [
        [0, 1, 2],
        [1, 2, 3],
        [2, 3, 4],
    ]
    y = [
        [0, 0, 0],
        [1, 1, 1],
        [2, 2, 2],
    ]
    grid = QuadModelGrid(x, y)
    assert grid.geometry is QuadGridGeometry.GENERAL


def test_from_general_but_raster():
    x = [
        [0, 1, 2],
        [0, 1, 2],
        [0, 1, 2],
    ]
    y = [
        [0, 0, 0],
        [1, 1, 1],
        [2, 2, 2],
    ]
    grid = QuadModelGrid(x, y)
    assert grid.geometry is QuadGridGeometry.RASTER


def test_from_general_but_rectilinear():
    x = [
        [0, 1, 2],
        [0, 1, 2],
        [0, 1, 2],
    ]
    y = [
        [0, 0, 0],
        [1, 1, 1],
        [3, 3, 3],
    ]
    grid = QuadModelGrid(x, y)
    assert grid.geometry is QuadGridGeometry.RECTILINEAR


def test_from_rectilinear():
    x = np.logspace(0, 10)
    y = np.logspace(10, 100)

    grid = QuadModelGrid.from_rectilinear(x, y)

    assert grid.shape == (50, 50)
    assert grid.geometry is QuadGridGeometry.RECTILINEAR


def test_xy_spacing_attribute():
    grid = QuadModelGrid.from_raster((4, 5), xy_spacing=(2.0, 1.0))
    assert grid.xy_spacing == (2.0, 1.0)


def test_xy_spacing_attribute_missing():
    x = np.logspace(0, 10)
    y = np.logspace(10, 100)

    grid = QuadModelGrid.from_rectilinear(x, y)
    with pytest.raises(AttributeError, match="xy_spacing is only defined for raster"):
        grid.xy_spacing


@pytest.mark.parametrize(
    "bc, expected",
    (
        ("open", NodeStatus.FIXED_VALUE),
        ("closed", NodeStatus.CLOSED),
        (NodeStatus.FIXED_VALUE, NodeStatus.FIXED_VALUE),
        (NodeStatus.CLOSED, NodeStatus.CLOSED),
    ),
)
def test_normalize_bc_at_edge_with_scalar(bc, expected):
    actual = _normalize_bc_at_edge(bc)
    assert actual == {
        "right": expected,
        "top": expected,
        "left": expected,
        "bottom": expected,
    }


def test_normalize_bc_at_edge_with_dict():
    actual = _normalize_bc_at_edge({"top": NodeStatus.CLOSED, "bottom": "closed"})
    assert actual == {
        "right": NodeStatus.FIXED_VALUE,
        "top": NodeStatus.CLOSED,
        "left": NodeStatus.FIXED_VALUE,
        "bottom": NodeStatus.CLOSED,
    }


@pytest.mark.parametrize("bc", (None, True, False, 0, 4))
def test_normalize_bc_at_edge_with_bad_bc(bc):
    with pytest.raises(TypeError, match="bc must be either a str"):
        _normalize_bc_at_edge(bc)


def test_make_node_status():
    actual = _make_node_status((3, 4), bc_at_edge={"top": 4})
    expected = [
        [0, 0, 0, 0],
        [0, 0, 0, 0],
        [4, 4, 4, 4],
    ]
    assert actual.dtype == np.uint8
    assert_array_equal(actual, np.asarray(expected))


@pytest.mark.parametrize("edge", ("TOP", "north", None))
def test_make_node_status_bad_edge(edge):
    with pytest.raises(KeyError, match="unknown value for edge"):
        _make_node_status((3, 4), bc_at_edge={edge: 4})

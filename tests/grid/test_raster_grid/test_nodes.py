import numpy as np
import pytest
from numpy.testing import assert_array_equal

from landlab import RasterModelGrid


def test_nodes_at_left_edge():
    """Test nodes at left edge of raster grid."""
    grid = RasterModelGrid((3, 4))
    assert_array_equal(grid.nodes_at_left_edge, np.array([0, 4, 8], dtype=np.int))

    vals_at_node = np.array([0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11])
    assert_array_equal(vals_at_node[grid.nodes_at_left_edge], np.array([0, 4, 8]))


def test_nodes_at_right_edge():
    """Test nodes at right edge of raster grid."""
    grid = RasterModelGrid((3, 4))
    assert_array_equal(grid.nodes_at_right_edge, np.array([3, 7, 11], dtype=np.int))


def test_nodes_at_top_edge():
    """Test nodes at top edge of raster grid."""
    grid = RasterModelGrid((3, 4))
    assert_array_equal(grid.nodes_at_top_edge, np.array([8, 9, 10, 11], dtype=np.int))


def test_nodes_at_bottom_edge():
    """Test nodes at bottom edge of raster grid."""
    grid = RasterModelGrid((3, 4))
    assert_array_equal(grid.nodes_at_bottom_edge, np.array([0, 1, 2, 3], dtype=np.int))


def test_nodes_at_edge():
    """Test nodes at edges (by name) of raster grid."""
    grid = RasterModelGrid((3, 4))
    for edge in ("right", "top", "left", "bottom"):
        assert_array_equal(
            grid.nodes_at_edge(edge), getattr(grid, "nodes_at_{0}_edge".format(edge))
        )
    with pytest.raises(ValueError):
        grid.nodes_at_edge("not-an-edge")


def test_grid_shape():
    """Test shape of grid."""
    grid = RasterModelGrid((3, 4))
    assert grid.shape == (3, 4)


def test_node_rows():
    """Test number of node rows."""
    grid = RasterModelGrid((4, 5))
    assert grid.number_of_node_rows == 4


def test_node_columns():
    """Test number of node columns."""
    grid = RasterModelGrid((4, 5))
    assert grid.number_of_node_columns == 5


def test_x_of_node():
    """Test x-coordinates of nodes."""
    grid = RasterModelGrid((3, 4), xy_spacing=(3.0, 2.0))
    assert_array_equal(
        grid.x_of_node,
        np.array([0.0, 3.0, 6.0, 9.0, 0.0, 3.0, 6.0, 9.0, 0.0, 3.0, 6.0, 9.0]),
    )

    with pytest.raises(ValueError):
        grid.x_of_node[0, 0] = 1.0


def test_y_of_node():
    """Test y-coordinates of nodes."""
    grid = RasterModelGrid((3, 4), xy_spacing=(3.0, 2.0))
    assert_array_equal(
        grid.y_of_node,
        np.array([0.0, 0.0, 0.0, 0.0, 2.0, 2.0, 2.0, 2.0, 4.0, 4.0, 4.0, 4.0]),
    )

    with pytest.raises(ValueError):
        grid.y_of_node[0, 0] = 1.0


def test_xy_of_node():
    """Test coordinates of nodes as x-y pairs."""
    grid = RasterModelGrid((3, 4), xy_spacing=(3.0, 2.0))
    assert_array_equal(
        grid.xy_of_node,
        np.array(
            [
                [0.0, 0.0],
                [3.0, 0.0],
                [6.0, 0.0],
                [9.0, 0.0],
                [0.0, 2.0],
                [3.0, 2.0],
                [6.0, 2.0],
                [9.0, 2.0],
                [0.0, 4.0],
                [3.0, 4.0],
                [6.0, 4.0],
                [9.0, 4.0],
            ]
        ),
    )
    assert_array_equal(grid.xy_of_node[:, 0], grid.x_of_node)
    assert_array_equal(grid.xy_of_node[:, 1], grid.y_of_node)

    with pytest.raises(ValueError):
        grid.xy_of_node[0, 0] = 1.0


def test_dx():
    """Test spacing of columns."""
    grid = RasterModelGrid((4, 5))
    assert grid.dx == 1.0

    grid = RasterModelGrid((4, 5), xy_spacing=2.0)
    assert grid.dx == 2.0

    grid = RasterModelGrid((4, 5), xy_spacing=(1.0, 2.0))
    assert grid.dx == 1.0


def test_dy():
    """Test spacing of rows."""
    grid = RasterModelGrid((4, 5))
    assert grid.dy == 1.0

    grid = RasterModelGrid((4, 5), xy_spacing=2.0)
    assert grid.dy == 2.0

    grid = RasterModelGrid((4, 5), xy_spacing=(1.0, 2.0))
    assert grid.dy == 2.0


def test_nodes_at_patch():
    """Test nodes around patches."""
    grid = RasterModelGrid((3, 3))
    assert_array_equal(
        grid.nodes_at_patch,
        np.array([[4, 3, 0, 1], [5, 4, 1, 2], [7, 6, 3, 4], [8, 7, 4, 5]], dtype=int),
    )
    with pytest.raises(ValueError):
        grid.nodes_at_patch[0, 0] = 42

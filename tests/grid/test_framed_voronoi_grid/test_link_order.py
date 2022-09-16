import numpy as np
from numpy.testing import assert_array_equal

from landlab import FramedVoronoiGrid


def test_rect_grid_link_order():
    """Test the order of links."""
    grid = FramedVoronoiGrid(
        (3, 2),
    )
    assert_array_equal(
        grid.nodes_at_link,
        [[0, 1], [0, 2], [2, 1], [1, 3], [2, 3], [2, 4], [4, 3], [3, 5], [4, 5]],
    )


def test_nodes_at_link():
    """Test nodes_at_link shares data with tail and head."""
    grid = FramedVoronoiGrid((3, 2))

    assert_array_equal(grid.nodes_at_link[:, 0], grid.node_at_link_tail)
    assert_array_equal(grid.nodes_at_link[:, 1], grid.node_at_link_head)

    assert np.may_share_memory(grid.nodes_at_link, grid.node_at_link_tail)
    assert np.may_share_memory(grid.nodes_at_link, grid.node_at_link_head)


def test_rect_face_at_link():
    grid = FramedVoronoiGrid((3, 3))
    assert_array_equal(
        grid.face_at_link,
        np.array([-1, -1, -1, 0, 1, 2, -1, 3, 4, -1, 5, 6, 7, -1, -1, -1]),
    )


def test_rect_link_at_face():
    grid = FramedVoronoiGrid((3, 3))
    assert_array_equal(grid.link_at_face, np.array([3, 4, 5, 7, 8, 10, 11, 12]))

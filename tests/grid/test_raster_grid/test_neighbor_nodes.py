import numpy as np
import pytest
from numpy.testing import assert_array_equal

from landlab import RasterModelGrid

X = RasterModelGrid.BAD_INDEX


def test_all_active_neighbors():
    rmg = RasterModelGrid((5, 4))
    expected = np.array(
        [
            [X, X, X, X],
            [X, 5, X, X],
            [X, 6, X, X],
            [X, X, X, X],
            [5, X, X, X],
            [6, 9, 4, 1],
            [7, 10, 5, 2],
            [X, X, 6, X],
            [9, X, X, X],
            [10, 13, 8, 5],
            [11, 14, 9, 6],
            [X, X, 10, X],
            [13, X, X, X],
            [14, 17, 12, 9],
            [15, 18, 13, 10],
            [X, X, 14, X],
            [X, X, X, X],
            [X, X, X, 13],
            [X, X, X, 14],
            [X, X, X, X],
        ]
    )
    assert_array_equal(rmg.active_adjacent_nodes_at_node, expected)


def test_all_neighbors():
    rmg = RasterModelGrid((5, 4))
    expected = np.array(
        [
            [1, 4, X, X],
            [2, 5, 0, X],
            [3, 6, 1, X],
            [X, 7, 2, X],
            [5, 8, X, 0],
            [6, 9, 4, 1],
            [7, 10, 5, 2],
            [X, 11, 6, 3],
            [9, 12, X, 4],
            [10, 13, 8, 5],
            [11, 14, 9, 6],
            [X, 15, 10, 7],
            [13, 16, X, 8],
            [14, 17, 12, 9],
            [15, 18, 13, 10],
            [X, 19, 14, 11],
            [17, X, X, 12],
            [18, X, 16, 13],
            [19, X, 17, 14],
            [X, X, 18, 15],
        ]
    )
    assert_array_equal(rmg.adjacent_nodes_at_node, expected)


def test_active_neighbor_list_with_scalar_arg():
    rmg = RasterModelGrid((5, 4))

    assert_array_equal(rmg.active_adjacent_nodes_at_node[6], np.array([7, 10, 5, 2]))
    assert_array_equal(rmg.active_adjacent_nodes_at_node[-1], np.array([X, X, X, X]))
    assert_array_equal(rmg.active_adjacent_nodes_at_node[-2], np.array([X, X, X, 14]))


def test_neighbor_list_with_scalar_arg():
    rmg = RasterModelGrid((5, 4))

    assert_array_equal(rmg.adjacent_nodes_at_node[6], np.array([7, 10, 5, 2]))
    assert_array_equal(rmg.adjacent_nodes_at_node[-1], np.array([X, X, 18, 15]))
    assert_array_equal(rmg.adjacent_nodes_at_node[-2], np.array([19, X, 17, 14]))


def test_active_neighbor_list_with_array_arg():
    rmg = RasterModelGrid((5, 4))
    assert_array_equal(
        rmg.active_adjacent_nodes_at_node[[6, -1]],
        np.array([[7, 10, 5, 2], [X, X, X, X]]),
    )


def test_neighbor_list_with_array_arg():
    rmg = RasterModelGrid((5, 4))
    assert_array_equal(
        rmg.adjacent_nodes_at_node[(6, -1), :],
        np.array([[7, 10, 5, 2], [X, X, 18, 15]]),
    )


def test_neighbor_list_is_read_only():
    rmg = RasterModelGrid((5, 4))
    with pytest.raises(ValueError):
        rmg.adjacent_nodes_at_node[0] = [1, 2, 3, 4]


def test_neighbors_is_contiguous():
    rmg = RasterModelGrid((5, 4))
    assert rmg.adjacent_nodes_at_node.flags["C_CONTIGUOUS"]


def test_active_neighbor_list_boundary():
    """All of the neighbor IDs for a boundary cell are -1."""
    rmg = RasterModelGrid((5, 4))
    import landlab.utils.structured_grid as sgrid

    rmg.status_at_node[(0, 1, 2, 3, 4, 7, 8, 11, 12, 15, 16, 17, 18, 19),] = (
        rmg.BC_NODE_IS_CLOSED
    )

    for node_id in sgrid.perimeter_iter(rmg.shape):
        assert_array_equal(
            rmg.active_adjacent_nodes_at_node[node_id], np.array([X, X, X, X])
        )


def test_all_diagonals():
    rmg = RasterModelGrid((5, 4))
    expected = np.array(
        [
            [5, X, X, X],
            [6, 4, X, X],
            [7, 5, X, X],
            [X, 6, X, X],
            [9, X, X, 1],
            [10, 8, 0, 2],
            [11, 9, 1, 3],
            [X, 10, 2, X],
            [13, X, X, 5],
            [14, 12, 4, 6],
            [15, 13, 5, 7],
            [X, 14, 6, X],
            [17, X, X, 9],
            [18, 16, 8, 10],
            [19, 17, 9, 11],
            [X, 18, 10, X],
            [X, X, X, 13],
            [X, X, 12, 14],
            [X, X, 13, 15],
            [X, X, 14, X],
        ]
    )
    assert_array_equal(rmg.diagonal_adjacent_nodes_at_node, expected)


def test_diagonal_list_with_scalar_arg():
    rmg = RasterModelGrid((5, 4))

    assert_array_equal(rmg.diagonal_adjacent_nodes_at_node[6], np.array([11, 9, 1, 3]))
    assert_array_equal(rmg.diagonal_adjacent_nodes_at_node[-1], np.array([X, X, 14, X]))
    assert_array_equal(
        rmg.diagonal_adjacent_nodes_at_node[-2], np.array([X, X, 13, 15])
    )


def test_diagonal_list_with_array_arg():
    rmg = RasterModelGrid((5, 4))
    assert_array_equal(
        rmg.diagonal_adjacent_nodes_at_node[(6, -1), :],
        np.array([[11, 9, 1, 3], [X, X, 14, X]]),
    )


def test_diagonal_list_is_read_only():
    rmg = RasterModelGrid((5, 4))
    with pytest.raises(ValueError):
        rmg.diagonal_adjacent_nodes_at_node[0] = [1, 2, 3, 4]


def test_diagonals_is_contiguous():
    rmg = RasterModelGrid((5, 4))
    assert rmg.diagonal_adjacent_nodes_at_node.flags["C_CONTIGUOUS"]

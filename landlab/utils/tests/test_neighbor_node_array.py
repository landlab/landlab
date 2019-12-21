import numpy as np
from numpy.testing import assert_array_equal

import landlab.utils.structured_grid as sgrid


def test_default():
    neighbors = sgrid.neighbor_node_array((2, 3))

    X = sgrid.BAD_INDEX_VALUE
    assert_array_equal(
        neighbors,
        np.array(
            [
                [1, 3, X, X],
                [2, 4, 0, X],
                [X, 5, 1, X],
                [4, X, X, 0],
                [5, X, 3, 1],
                [X, X, 4, 2],
            ]
        ).T,
    )

    assert neighbors.flags["C_CONTIGUOUS"]
    assert neighbors.base is None


def test_set_out_of_bounds():
    neighbors = sgrid.neighbor_node_array((2, 3), inactive=-1)
    assert_array_equal(
        neighbors,
        np.array(
            [
                [1, 3, -1, -1],
                [2, 4, 0, -1],
                [-1, 5, 1, -1],
                [4, -1, -1, 0],
                [5, -1, 3, 1],
                [-1, -1, 4, 2],
            ]
        ).T,
    )


def test_boundary_node_mask_no_actives():
    neighbors = sgrid.neighbor_node_array(
        (2, 3), inactive=-2, closed_boundary_nodes=np.arange(6)
    )
    assert_array_equal(neighbors, -2 * np.ones((6, 4)).T)


def test_boundary_node_mask_3x3():
    neighbors = sgrid.neighbor_node_array(
        (3, 3), inactive=-2, closed_boundary_nodes=[0, 1, 2, 3, 5, 6, 7, 8]
    )
    assert_array_equal(neighbors, -2 * np.ones((9, 4)).T)


def test_boundary_node_mask():
    X = -1
    neighbors = sgrid.neighbor_node_array(
        (3, 3), inactive=X, closed_boundary_nodes=[0, 2, 3, 5, 6, 7, 8]
    )
    expected = np.array(
        [
            [X, X, X, X],
            [X, 4, X, X],
            [X, X, X, X],
            [X, X, X, X],
            [X, X, X, 1],
            [X, X, X, X],
            [X, X, X, X],
            [X, X, X, X],
            [X, X, X, X],
        ]
    )
    assert_array_equal(neighbors, expected.T)


def test_open_boundary():
    X = -1
    neighbors = sgrid.neighbor_node_array(
        (3, 3),
        closed_boundary_nodes=[3, 5, 6, 7, 8],
        open_boundary_nodes=[0, 1, 2],
        inactive=X,
    )
    expected = np.array(
        [
            [X, X, X, X],
            [X, 4, X, X],
            [X, X, X, X],
            [X, X, X, X],
            [X, X, X, 1],
            [X, X, X, X],
            [X, X, X, X],
            [X, X, X, X],
            [X, X, X, X],
        ]
    )
    assert_array_equal(neighbors.T, expected)


def test_open_boundary_4x5():
    X = -1
    neighbors = sgrid.neighbor_node_array(
        (4, 5),
        closed_boundary_nodes=[5, 10, 15, 16, 17, 18, 19],
        open_boundary_nodes=[0, 1, 2, 3, 4, 9, 14],
        inactive=X,
    )
    expected = np.array(
        [
            [X, X, X, X],
            [X, 6, X, X],
            [X, 7, X, X],
            [X, 8, X, X],
            [X, X, X, X],
            [X, X, X, X],
            [7, 11, X, 1],
            [8, 12, 6, 2],
            [9, 13, 7, 3],
            [X, X, 8, X],
            [X, X, X, X],
            [12, X, X, 6],
            [13, X, 11, 7],
            [14, X, 12, 8],
            [X, X, 13, X],
            [X, X, X, X],
            [X, X, X, X],
            [X, X, X, X],
            [X, X, X, X],
            [X, X, X, X],
        ]
    )
    assert_array_equal(neighbors.T, expected)


def test_open_boundary_on_perimeter():
    X = -1

    neighbors = sgrid.neighbor_node_array(
        (5, 4), open_boundary_nodes=sgrid.boundary_nodes((5, 4)), inactive=X
    )

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
    assert_array_equal(neighbors, expected.T)


def test_closed_boundary_on_perimeter():
    X = -1

    neighbors = sgrid.neighbor_node_array(
        (5, 4), closed_boundary_nodes=sgrid.boundary_nodes((5, 4)), inactive=X
    )

    expected = np.array(
        [
            [X, X, X, X],
            [X, X, X, X],
            [X, X, X, X],
            [X, X, X, X],
            [X, X, X, X],
            [6, 9, X, X],
            [X, 10, 5, X],
            [X, X, X, X],
            [X, X, X, X],
            [10, 13, X, 5],
            [X, 14, 9, 6],
            [X, X, X, X],
            [X, X, X, X],
            [14, X, X, 9],
            [X, X, 13, 10],
            [X, X, X, X],
            [X, X, X, X],
            [X, X, X, X],
            [X, X, X, X],
            [X, X, X, X],
        ]
    )
    assert_array_equal(neighbors, expected.T)


def test_closed_boundaries_in_middle():
    X = -1

    neighbors = sgrid.neighbor_node_array(
        (5, 4), closed_boundary_nodes=sgrid.boundary_nodes((5, 4)), inactive=X
    )

    expected = np.array(
        [
            [X, X, X, X],
            [X, X, X, X],
            [X, X, X, X],
            [X, X, X, X],
            [X, X, X, X],
            [6, 9, X, X],
            [X, 10, 5, X],
            [X, X, X, X],
            [X, X, X, X],
            [10, 13, X, 5],
            [X, 14, 9, 6],
            [X, X, X, X],
            [X, X, X, X],
            [14, X, X, 9],
            [X, X, 13, 10],
            [X, X, X, X],
            [X, X, X, X],
            [X, X, X, X],
            [X, X, X, X],
            [X, X, X, X],
        ]
    )
    assert_array_equal(neighbors, expected.T)

#! /usr/bin/env python
import numpy as np

from nose.tools import assert_is, assert_is_instance
from numpy.testing import assert_array_equal

from landlab.grid.diagonals import create_nodes_at_diagonal


def test_nodes_at_diagonal():
    """Test tail and head nodes of diagonals."""
    diagonals = create_nodes_at_diagonal((2, 2))
    assert_array_equal(diagonals, [[0, 3], [1, 2]])

    diagonals = create_nodes_at_diagonal((4, 3))
    assert_is_instance(diagonals, np.ndarray)
    assert_array_equal(diagonals, [[0,  4], [1, 3], [1,  5], [2,  4],
                                   [3,  7], [4, 6], [4,  8], [5,  7],
                                   [6, 10], [7, 9], [7, 11], [8, 10]])


def test_nodes_at_diagonal_1d():
    """Test nodes at diagonals for 1d grid."""
    diagonals = create_nodes_at_diagonal((1, 2))
    assert_array_equal(diagonals, np.array([], dtype=int).reshape((0, 2)))

    diagonals = create_nodes_at_diagonal((4, 1))
    assert_array_equal(diagonals, np.array([], dtype=int).reshape((0, 2)))


def test_nodes_at_diagonal_out_keyword():
    """Test out keyword for nodes_at_diagonal."""
    buffer = np.empty((4, 2), dtype=int)
    diagonals = create_nodes_at_diagonal((3, 2), out=buffer)

    assert_is(buffer, diagonals)
    assert_array_equal(diagonals, [[0, 3], [1, 2], [2, 5], [3, 4]])

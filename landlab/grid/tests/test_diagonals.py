#! /usr/bin/env python
import numpy as np

from nose.tools import (assert_is, assert_is_instance, assert_tuple_equal,
                        assert_raises, assert_false)
from numpy.testing import assert_array_equal

from landlab import RasterModelGrid
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


def test_values_are_cached():
    """Test that attributes of diagonals are cached."""
    names = (
        'diagonals_at_node',
        'diagonal_dirs_at_node',
        'diagonal_adjacent_nodes_at_node',
        'd8_adjacent_nodes_at_node',
        'nodes_at_diagonal',
        'nodes_at_d8',
        'd8s_at_node',
        'd8_dirs_at_node',
        # 'd8_status_at_node',
        'length_of_diagonal',
        'length_of_d8',
        'status_at_diagonal',
        'diagonal_status_at_node',
        'active_diagonals',
        'active_diagonal_dirs_at_node',
        'status_at_d8',
        'active_d8',
        'active_d8_dirs_at_node',
    )

    for name in names:
        def _check_value_is_cached(attr):
            grid = RasterModelGrid((4, 3))
            x = getattr(grid, attr)
            assert_is(x, getattr(grid, attr))
        _check_value_is_cached.description = 'Test {name} is cached'.format(
            name=name)
        yield _check_value_is_cached, name


def test_values_are_readonly():
    """Test that diagonals attributes are readonly."""
    names = (
        'diagonals_at_node',
        'diagonal_dirs_at_node',
        'diagonal_adjacent_nodes_at_node',
        'd8_adjacent_nodes_at_node',
        'nodes_at_diagonal',
        'nodes_at_d8',
        'd8s_at_node',
        'd8_dirs_at_node',
        # 'd8_status_at_node',
        'length_of_diagonal',
        'length_of_d8',
        'status_at_diagonal',
        'diagonal_status_at_node',
        'active_diagonals',
        'active_diagonal_dirs_at_node',
        'status_at_d8',
        'active_d8',
        'active_d8_dirs_at_node',
    )

    for name in names:
        def _check_value_is_readonly(attr):
            grid = RasterModelGrid((4, 3))
            x = getattr(grid, attr)
            with assert_raises(ValueError):
                x[0] = 999
            assert_false(x.flags['WRITEABLE'])
        _check_value_is_readonly.description = 'Test {name} is readonly'.format(
            name=name)
        yield _check_value_is_readonly, name


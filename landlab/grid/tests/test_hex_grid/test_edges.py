import numpy as np
from nose.tools import assert_false, assert_is, assert_raises
from numpy.testing import assert_array_equal

from landlab import HexModelGrid


def test_perimeter_nodes():
    """Test perimeter nodes of a hex grid."""
    grid = HexModelGrid(3, 4, shape="rect")
    assert_array_equal(grid.perimeter_nodes, [3, 7, 11, 10, 9, 8, 4, 0, 1, 2])


def test_right_edge_nodes():
    """Test right edge nodes of a hex grid."""
    grid = HexModelGrid(3, 4, shape="rect")
    assert_array_equal(grid.nodes_at_right_edge, [3, 7, 11])


def test_top_edge_nodes():
    """Test top edge nodes of a hex grid."""
    grid = HexModelGrid(3, 4, shape="rect")
    assert_array_equal(grid.nodes_at_top_edge, [8, 9, 10, 11])


def test_left_edge_nodes():
    """Test left edge nodes of a hex grid."""
    grid = HexModelGrid(3, 4, shape="rect")
    assert_array_equal(grid.nodes_at_left_edge, [0, 4, 8])


def test_bottom_edge_nodes():
    """Test bottom edge nodes of a hex grid."""
    grid = HexModelGrid(3, 4, shape="rect")
    assert_array_equal(grid.nodes_at_bottom_edge, [0, 1, 2, 3])


def test_edges_are_readonly():
    names = [
        "nodes_at_right_edge",
        "nodes_at_top_edge",
        "nodes_at_left_edge",
        "nodes_at_bottom_edge",
    ]

    for name in names:

        def _test_readonly():
            grid = HexModelGrid(3, 4, shape="rect")
            assert_false(grid.perimeter_nodes.flags["WRITEABLE"])
            with assert_raises(ValueError):
                getattr(grid, name)[0] = 999

        _test_readonly.description = "Test {name} is readonly".format(name=name)
        yield _test_readonly


def test_edges_are_cached():
    names = [
        "nodes_at_right_edge",
        "nodes_at_top_edge",
        "nodes_at_left_edge",
        "nodes_at_bottom_edge",
    ]

    for name in names:

        def _test_cached():
            grid = HexModelGrid(3, 4, shape="rect")
            x = grid.perimeter_nodes
            assert_is(grid.perimeter_nodes, x)

        _test_cached.description = "Test {name} is cached".format(name=name)
        yield _test_cached

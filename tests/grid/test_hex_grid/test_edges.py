import numpy as np
import pytest

from landlab import HexModelGrid


def test_perimeter_nodes():
    """Test perimeter nodes of a hex grid."""
    grid = HexModelGrid((3, 4), node_layout="rect")
    assert np.all(grid.perimeter_nodes == [3, 7, 11, 10, 9, 8, 4, 0, 1, 2])


def test_right_edge_nodes():
    """Test right edge nodes of a hex grid."""
    grid = HexModelGrid((3, 4), node_layout="rect")
    assert np.all(grid.nodes_at_right_edge == [3, 7, 11])


def test_top_edge_nodes():
    """Test top edge nodes of a hex grid."""
    grid = HexModelGrid((3, 4), node_layout="rect")
    assert np.all(grid.nodes_at_top_edge == [8, 9, 10, 11])


def test_left_edge_nodes():
    """Test left edge nodes of a hex grid."""
    grid = HexModelGrid((3, 4), node_layout="rect")
    assert np.all(grid.nodes_at_left_edge == [0, 4, 8])


def test_bottom_edge_nodes():
    """Test bottom edge nodes of a hex grid."""
    grid = HexModelGrid((3, 4), node_layout="rect")
    assert np.all(grid.nodes_at_bottom_edge == [0, 1, 2, 3])


def test_edges_are_readonly(edge_name):
    grid = HexModelGrid((3, 4), node_layout="rect")
    assert not grid.perimeter_nodes.flags["WRITEABLE"]
    with pytest.raises(ValueError):
        getattr(grid, "nodes_at_" + edge_name)[0] = 999


def test_edges_are_cached(edge_name):
    grid = HexModelGrid((3, 4), node_layout="rect")
    x = grid.perimeter_nodes
    assert grid.perimeter_nodes is x
    x = getattr(grid, "nodes_at_" + edge_name)
    assert getattr(grid, "nodes_at_" + edge_name) is x

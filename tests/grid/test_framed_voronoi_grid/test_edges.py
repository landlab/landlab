import numpy as np
import pytest

from landlab import FramedVoronoiGrid


def test_rect_perimeter_nodes():
    """Test perimeter nodes of a framed voronoi grid with rectangular layout."""
    grid = FramedVoronoiGrid((3, 4))
    assert np.all(grid.perimeter_nodes == [3, 7, 11, 10, 9, 8, 4, 0, 1, 2])


def test_rect_right_edge_nodes():
    """Test right edge nodes of a framed voronoi grid with rectangular layout."""
    grid = FramedVoronoiGrid((3, 4))
    assert np.all(grid.nodes_at_right_edge == [3, 7, 11])


def test_rect_top_edge_nodes():
    """Test top edge nodes of a framed voronoi grid with rectangular layout."""
    grid = FramedVoronoiGrid((3, 4))
    assert np.all(grid.nodes_at_top_edge == [8, 9, 10, 11])


def test_rect_left_edge_nodes():
    """Test left edge nodes of a framed voronoi grid with rectangular layout."""
    grid = FramedVoronoiGrid((3, 4))
    assert np.all(grid.nodes_at_left_edge == [0, 4, 8])


def test_rect_bottom_edge_nodes():
    """Test bottom edge nodes of framed voronoi grid with rectangular layout."""
    grid = FramedVoronoiGrid((3, 4))
    assert np.all(grid.nodes_at_bottom_edge == [0, 1, 2, 3])


def test_edges_are_readonly():
    grid = FramedVoronoiGrid((3, 4))
    edge_name = "bottom_edge"
    assert not grid.perimeter_nodes.flags["WRITEABLE"]
    with pytest.raises(ValueError):
        getattr(grid, "nodes_at_" + edge_name)[0] = 999

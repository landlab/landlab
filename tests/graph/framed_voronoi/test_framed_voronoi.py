import numpy as np
import pytest
from hypothesis import given
from hypothesis.strategies import integers
from hypothesis.strategies import lists
from numpy.testing import assert_array_almost_equal
from numpy.testing import assert_array_equal
from numpy.testing import assert_array_less
from pytest import approx

from landlab.graph.framed_voronoi.framed_voronoi import FramedVoronoiGraph
from landlab.graph.framed_voronoi.framed_voronoi import HorizontalRectVoronoiGraph


def test_number_of_nodes_horizontal_rect():
    assert HorizontalRectVoronoiGraph.number_of_nodes((1, 2)) == 2
    assert HorizontalRectVoronoiGraph.number_of_nodes((1, 3)) == 3
    assert HorizontalRectVoronoiGraph.number_of_nodes((2, 2)) == 4
    assert HorizontalRectVoronoiGraph.number_of_nodes((2, 3)) == 6
    assert HorizontalRectVoronoiGraph.number_of_nodes((3, 2)) == 6
    assert HorizontalRectVoronoiGraph.number_of_nodes((3, 3)) == 9


@pytest.mark.parametrize("n_rows", (3,))
@pytest.mark.parametrize("at", ("nodes", "links", "patches"))
def test_create_rect_graph(n_rows, at):
    expected = {"nodes": 6, "links": 9, "patches": 4}
    shape = (n_rows, 2)
    graph = FramedVoronoiGraph(shape, sort=True)
    assert getattr(graph, f"number_of_{at}") == expected[at]


def test_create_rect():
    """Test creating a hex graph with rectangular layout."""
    graph = FramedVoronoiGraph((3, 2), sort=True)

    assert graph.number_of_nodes == 6
    assert graph.number_of_links == 9
    assert graph.number_of_patches == 4


@given(shape=lists(integers(min_value=3, max_value=5), min_size=2, max_size=2))  # 32
def test_spacing(shape):
    """Test spacing of nodes."""
    graph = FramedVoronoiGraph(shape, xy_spacing=1.0, xy_min_spacing=0.5)
    assert_array_less(0.497, graph.length_of_link)
    assert_array_less(graph.length_of_link, np.sqrt(1.5**2 + 1.5**2))


@given(shape=lists(integers(min_value=3, max_value=32), min_size=2, max_size=2))
def test_origin_keyword(shape):
    """Test setting the origin."""
    graph = FramedVoronoiGraph(shape)

    assert np.min(graph.x_of_node) == approx(0.0)
    assert np.min(graph.y_of_node) == approx(0.0)

    graph = FramedVoronoiGraph(shape, xy_of_lower_left=(0.5, 0.25))

    assert np.min(graph.x_of_node[0]) == approx(0.5)
    assert np.min(graph.y_of_node[0]) == approx(0.25)


def test_rect_orientation():
    """Test horizontal orientation."""
    graph = FramedVoronoiGraph(
        (3, 3),
    )
    assert_array_almost_equal(
        graph.x_of_node[0:4],
        [0.0, 1.0, 2.0, 0.0],
    )
    assert_array_almost_equal(
        graph.x_of_node[5:9],
        [2.0, 0.0, 1.0, 2.0],
    )


def test_perimeter_nodes_rect():
    graph = FramedVoronoiGraph(
        (3, 4),
    )
    assert_array_equal(graph.perimeter_nodes, [3, 7, 11, 10, 9, 8, 4, 0, 1, 2])


def test_rect_adjacent_nodes_at_node():
    graph = FramedVoronoiGraph(
        (3, 3),
        sort=True,
    )
    assert_array_equal(
        graph.adjacent_nodes_at_node,
        [
            [1, 3, -1, -1, -1, -1],
            [2, 4, 3, 0, -1, -1],
            [5, 4, 1, -1, -1, -1],
            [4, 6, 0, 1, -1, -1],
            [5, 7, 6, 3, 1, 2],
            [8, 7, 4, 2, -1, -1],
            [7, 3, 4, -1, -1, -1],
            [8, 6, 4, 5, -1, -1],
            [7, 5, -1, -1, -1, -1],
        ],
    )


def test_rect_patches_at_node():
    graph = FramedVoronoiGraph(
        (3, 3),
        sort=True,
    )
    assert np.array_equal(
        graph.patches_at_node,
        [
            [2, -1, -1, -1, -1, -1],
            [1, 3, 2, -1, -1, -1],
            [5, 1, -1, -1, -1, -1],
            [0, 2, 3, -1, -1, -1],
            [7, 4, 0, 3, 1, 5],
            [6, 7, 5, -1, -1, -1],
            [0, 4, -1, -1, -1, -1],
            [4, 7, 6, -1, -1, -1],
            [6, -1, -1, -1, -1, -1],
        ],
    )


@pytest.mark.parametrize("n_cols", [(3)])
@pytest.mark.parametrize("n_rows", (1, 3))
def test_xy_of_node_rect_horizontal(n_rows, n_cols):
    expected = {
        (1, 3): ([0.0, 1.0, 2.0], [0.0, 0.0, 0.0]),
        (3, 3): (
            [0.0, 1.0, 2.0],
            [0.0, 0.0, 0.0],
        ),
    }
    x_of_node, y_of_node = HorizontalRectVoronoiGraph.xy_of_node((n_rows, n_cols))

    assert np.all(x_of_node[0:3] == approx(expected[(n_rows, n_cols)][0]))
    assert np.all(y_of_node[0:3] == approx(expected[(n_rows, n_cols)][1]))


@pytest.mark.parametrize("n_cols", (2, 3))
@pytest.mark.parametrize("n_rows", (1, 2, 3))
def test_xy_of_node_lower_left(layout_graph, n_rows, n_cols):
    (x_of_node, y_of_node) = layout_graph.xy_of_node((n_rows, n_cols))

    assert np.min(x_of_node) == approx(0.0)
    assert np.min(y_of_node) == approx(0.0)

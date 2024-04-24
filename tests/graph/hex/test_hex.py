import numpy as np
import pytest
from hypothesis import given
from hypothesis.strategies import integers
from hypothesis.strategies import lists
from numpy.testing import assert_array_almost_equal
from numpy.testing import assert_array_equal
from pytest import approx

from landlab.graph import TriGraph
from landlab.graph.hex.hex import HorizontalHexTriGraph
from landlab.graph.hex.hex import HorizontalRectTriGraph
from landlab.graph.hex.hex import VerticalHexTriGraph
from landlab.graph.hex.hex import VerticalRectTriGraph


def test_number_of_nodes_horizontal_rect():
    assert HorizontalRectTriGraph.number_of_nodes((1, 2)) == 2
    assert HorizontalRectTriGraph.number_of_nodes((1, 3)) == 3
    assert HorizontalRectTriGraph.number_of_nodes((2, 2)) == 4
    assert HorizontalRectTriGraph.number_of_nodes((2, 3)) == 6
    assert HorizontalRectTriGraph.number_of_nodes((3, 2)) == 6
    assert HorizontalRectTriGraph.number_of_nodes((3, 3)) == 9


def test_number_of_nodes_vertical_rect():
    assert VerticalRectTriGraph.number_of_nodes((1, 2)) == 2
    assert VerticalRectTriGraph.number_of_nodes((1, 3)) == 3
    assert VerticalRectTriGraph.number_of_nodes((2, 2)) == 4
    assert VerticalRectTriGraph.number_of_nodes((2, 3)) == 6
    assert VerticalRectTriGraph.number_of_nodes((3, 2)) == 6
    assert VerticalRectTriGraph.number_of_nodes((3, 3)) == 9


def test_number_of_nodes_horizontal_hex():
    assert HorizontalHexTriGraph.number_of_nodes((1, 2)) == 2
    assert HorizontalHexTriGraph.number_of_nodes((1, 3)) == 3
    assert HorizontalHexTriGraph.number_of_nodes((2, 2)) == 5
    assert HorizontalHexTriGraph.number_of_nodes((2, 3)) == 7
    assert HorizontalHexTriGraph.number_of_nodes((3, 2)) == 7
    assert HorizontalHexTriGraph.number_of_nodes((3, 3)) == 10


def test_number_of_nodes_vertical_hex():
    assert VerticalHexTriGraph.number_of_nodes((1, 2)) == 3
    assert VerticalHexTriGraph.number_of_nodes((1, 3)) == 4
    assert VerticalHexTriGraph.number_of_nodes((2, 2)) == 5
    assert VerticalHexTriGraph.number_of_nodes((2, 3)) == 7
    assert VerticalHexTriGraph.number_of_nodes((3, 2)) == 7
    assert VerticalHexTriGraph.number_of_nodes((3, 3)) == 10


@given(shape=lists(integers(min_value=3, max_value=1024), min_size=2, max_size=2))
def test_number_of_nodes_symetric_rect(shape):
    assert VerticalRectTriGraph.number_of_nodes(
        shape
    ) == HorizontalRectTriGraph.number_of_nodes(shape[::-1])


@given(shape=lists(integers(min_value=3, max_value=1024), min_size=2, max_size=2))
def test_number_of_nodes_symetric_hex(shape):
    assert VerticalHexTriGraph.number_of_nodes(
        shape
    ) == HorizontalHexTriGraph.number_of_nodes(shape[::-1])


@pytest.mark.parametrize("n_rows", (3,))
@pytest.mark.parametrize("node_layout", ("rect", "hex"))
@pytest.mark.parametrize("orientation", ("horizontal", "vertical"))
@pytest.mark.parametrize("at", ("nodes", "links", "patches"))
def test_create_hex_graph(n_rows, node_layout, orientation, at):
    expected = {
        "rect": {
            "horizontal": {"nodes": 6, "links": 9, "patches": 4},
            "vertical": {"nodes": 6, "links": 9, "patches": 4},
        },
        "hex": {
            "horizontal": {"nodes": 7, "links": 12, "patches": 6},
            "vertical": {"nodes": 7, "links": 12, "patches": 6},
        },
    }
    if orientation == "vertical":
        shape = (2, n_rows)
    else:
        shape = (n_rows, 2)
    graph = TriGraph(shape, node_layout=node_layout, orientation=orientation, sort=True)
    assert getattr(graph, f"number_of_{at}") == expected[node_layout][orientation][at]


def test_create_rect():
    """Test creating a hex graph with rectangular layout."""
    graph = TriGraph((3, 2), node_layout="rect", sort=True)

    assert graph.number_of_nodes == 6
    assert graph.number_of_links == 9
    assert graph.number_of_patches == 4


def test_create_hex():
    """Test creating a hex graph with hex layout."""
    graph = TriGraph((3, 2), node_layout="hex", sort=True)

    assert graph.number_of_nodes == 7
    assert graph.number_of_links == 12
    assert graph.number_of_patches == 6


@given(shape=lists(integers(min_value=3, max_value=32), min_size=2, max_size=2))
def test_spacing(shape):
    """Test spacing of nodes."""
    graph = TriGraph(shape)
    assert_array_almost_equal(graph.length_of_link, 1.0)

    graph = TriGraph(shape, spacing=2)
    assert_array_almost_equal(graph.length_of_link, 2.0)


@given(shape=lists(integers(min_value=3, max_value=32), min_size=2, max_size=2))
@pytest.mark.parametrize("orientation", ("horizontal", "vertical"))
@pytest.mark.parametrize("node_layout", ("hex", "rect"))
def test_origin_keyword(node_layout, orientation, shape):
    """Test setting the origin."""
    graph = TriGraph(shape)

    assert np.min(graph.x_of_node) == approx(0.0)
    assert np.min(graph.y_of_node) == approx(0.0)

    graph = TriGraph(shape, xy_of_lower_left=(0.5, 0.25))

    assert np.min(graph.x_of_node[0]) == approx(0.5)
    assert np.min(graph.y_of_node[0]) == approx(0.25)


def test_orientation():
    """Test vertical and horizontal orientation."""
    graph = TriGraph((3, 3), orientation="vertical")
    assert_array_almost_equal(
        graph.y_of_node, [0.0, 0.0, 0.5, 1.0, 1.0, 1.5, 2.0, 2.0, 2.5]
    )

    graph = TriGraph((3, 3), orientation="horizontal")
    assert_array_almost_equal(
        graph.x_of_node, [0.0, 1.0, 2.0, 0.5, 1.5, 2.5, 0.0, 1.0, 2.0]
    )


def test_perimeter_nodes_rect():
    graph = TriGraph((3, 4), node_layout="rect")
    assert_array_equal(graph.perimeter_nodes, [3, 7, 11, 10, 9, 8, 4, 0, 1, 2])


def test_perimeter_nodes_hex():
    graph = TriGraph((4, 2), node_layout="hex")
    assert_array_equal(graph.perimeter_nodes, [8, 11, 10, 9, 5, 2, 0, 1, 4])


def test_adjacent_nodes_at_node():
    graph = TriGraph((3, 3), node_layout="hex", sort=True)
    assert_array_equal(
        graph.adjacent_nodes_at_node,
        [
            [1, 4, 3, -1, -1, -1],
            [2, 5, 4, 0, -1, -1],
            [6, 5, 1, -1, -1, -1],
            [4, 7, 0, -1, -1, -1],
            [5, 8, 7, 3, 0, 1],
            [6, 9, 8, 4, 1, 2],
            [9, 5, 2, -1, -1, -1],
            [8, 3, 4, -1, -1, -1],
            [9, 7, 4, 5, -1, -1],
            [8, 5, 6, -1, -1, -1],
        ],
    )


def test_patches_at_node():
    grid = TriGraph((3, 3), node_layout="hex", sort=True)
    assert_array_equal(
        grid.patches_at_node,
        [
            [0, 2, -1, -1, -1, -1],
            [1, 3, 0, -1, -1, -1],
            [4, 1, -1, -1, -1, -1],
            [5, 2, -1, -1, -1, -1],
            [6, 8, 5, 2, 0, 3],
            [7, 9, 6, 3, 1, 4],
            [7, 4, -1, -1, -1, -1],
            [5, 8, -1, -1, -1, -1],
            [8, 6, 9, -1, -1, -1],
            [9, 7, -1, -1, -1, -1],
        ],
    )


@pytest.mark.parametrize("n_cols", (2, 3))
@pytest.mark.parametrize("n_rows", (2, 3))
def test_xy_of_node_rect_vertical(n_rows, n_cols):
    expected = {
        (2, 2): ([0, 1, 0, 1], [0, 0.5, 1, 1.5]),
        (2, 3): ([0, 2, 1, 0, 2, 1], [0, 0, 0.5, 1, 1, 1.5]),
        (3, 2): ([0, 1, 0, 1, 0, 1], [0, 0.5, 1, 1.5, 2, 2.5]),
        (3, 3): ([0, 2, 1, 0, 2, 1, 0, 2, 1], [0, 0, 0.5, 1, 1, 1.5, 2, 2, 2.5]),
    }
    x_of_node, y_of_node = VerticalRectTriGraph.xy_of_node((n_rows, n_cols))

    assert np.all(
        x_of_node / np.sin(np.pi / 3.0) == approx(expected[(n_rows, n_cols)][0])
    )
    assert np.all(y_of_node == approx(expected[(n_rows, n_cols)][1]))


@pytest.mark.parametrize("n_cols", (2, 3))
@pytest.mark.parametrize("n_rows", (1, 2, 3))
def test_xy_of_node_hex_vertical(n_rows, n_cols):
    expected = {
        (1, 2): ([1.0, 0, 1.0], [0, 0.5, 1]),
        (1, 3): ([1.0, 0, 2, 1.0], [0, 0.5, 0.5, 1]),
        (2, 2): ([1.0, 0, 1.0, 0, 1.0], [0, 0.5, 1, 1.5, 2]),
        (2, 3): ([1.0, 0, 2, 1.0, 0, 2, 1.0], [0, 0.5, 0.5, 1, 1.5, 1.5, 2]),
        (3, 2): ([1.0, 0, 1.0, 0, 1.0, 0.0, 1.0], [0, 0.5, 1, 1.5, 2, 2.5, 3]),
        (3, 3): (
            [1.0, 0, 2, 1.0, 0, 2, 1.0, 0, 2, 1.0],
            [0, 0.5, 0.5, 1, 1.5, 1.5, 2, 2.5, 2.5, 3],
        ),
    }
    x_of_node, y_of_node = VerticalHexTriGraph.xy_of_node((n_rows, n_cols))

    assert np.all(
        x_of_node / np.sin(np.pi / 3.0) == approx(expected[(n_rows, n_cols)][0])
    )
    assert np.all(y_of_node == approx(expected[(n_rows, n_cols)][1]))


def test_xy_of_node_spacing(hex_layout):
    x_of_node_expected, y_of_node_expected = hex_layout.xy_of_node((3, 4))
    x_of_node, y_of_node = hex_layout.xy_of_node((3, 4), spacing=2.0)

    assert_array_almost_equal(x_of_node / 2.0, x_of_node_expected)
    assert_array_almost_equal(y_of_node / 2.0, y_of_node_expected)


@pytest.mark.parametrize("n_cols", (2, 3))
@pytest.mark.parametrize("n_rows", (1, 2, 3))
def test_xy_of_node_lower_left(hex_layout, n_rows, n_cols):
    (x_of_node, y_of_node) = hex_layout.xy_of_node((n_rows, n_cols))

    assert np.min(x_of_node) == approx(0.0)
    assert np.min(y_of_node) == approx(0.0)

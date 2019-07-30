import numpy as np
import pytest
from numpy.testing import assert_array_almost_equal, assert_array_equal
from pytest import approx

from landlab.graph import TriGraph
from landlab.graph.hex.ext.hex import (
    fill_xy_of_node_hex_vertical,
    fill_xy_of_node_rect_vertical,
)
from landlab.graph.hex.hex import number_of_nodes, setup_xy_of_node


@pytest.mark.parametrize("n_cols", (2, 3))
@pytest.mark.parametrize("n_rows", (1, 2, 3))
@pytest.mark.parametrize("orientation", ("horizontal", "vertical"))
@pytest.mark.parametrize("node_layout", ("rect", "hex"))
def test_number_of_nodes(node_layout, orientation, n_rows, n_cols):
    expected = {
        "rect": {
            "horizontal": {
                (1, 2): 2,
                (1, 3): 3,
                (2, 2): 4,
                (2, 3): 6,
                (3, 2): 6,
                (3, 3): 9,
            },
            "vertical": {
                (1, 2): 2,
                (1, 3): 3,
                (2, 2): 4,
                (2, 3): 6,
                (3, 2): 6,
                (3, 3): 9,
            },
        },
        "hex": {
            "horizontal": {
                (1, 2): 2,
                (1, 3): 3,
                (2, 2): 5,
                (2, 3): 7,
                (3, 2): 7,
                (3, 3): 10,
            },
            "vertical": {
                (1, 2): 3,
                (1, 3): 4,
                (2, 2): 5,
                (2, 3): 7,
                (3, 2): 7,
                (3, 3): 10,
            },
        },
    }
    assert (
        number_of_nodes(
            (n_rows, n_cols), orientation=orientation, node_layout=node_layout
        )
        == expected[node_layout][orientation][(n_rows, n_cols)]
    )


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
    graph = TriGraph(shape, node_layout=node_layout, orientation=orientation)
    assert (
        getattr(graph, "number_of_{at}".format(at=at))
        == expected[node_layout][orientation][at]
    )


def test_create_rect():
    """Test creating a hex graph with rectangular layout."""
    graph = TriGraph((3, 2), node_layout="rect")

    assert graph.number_of_nodes == 6
    assert graph.number_of_links == 9
    assert graph.number_of_patches == 4


def test_create_hex():
    """Test creating a hex graph with hex layout."""
    graph = TriGraph((3, 2), node_layout="hex")

    assert graph.number_of_nodes == 7
    assert graph.number_of_links == 12
    assert graph.number_of_patches == 6


@pytest.mark.skip
def test_create_rect1():
    """Test creating a hex graph."""
    graph = TriGraph((3, 2), node_layout="rect1")

    assert graph.number_of_nodes == 7
    assert graph.number_of_links == 12
    assert graph.number_of_patches == 6


def test_spacing():
    """Test spacing of nodes."""
    graph = TriGraph((20, 31))
    assert_array_almost_equal(graph.length_of_link, 1.0)

    graph = TriGraph((31, 20), spacing=2)
    assert_array_almost_equal(graph.length_of_link, 2.0)


@pytest.mark.parametrize("n_cols", (3, 4))
@pytest.mark.parametrize("n_rows", (3, 4))
@pytest.mark.parametrize("orientation", ("horizontal", "vertical"))
@pytest.mark.parametrize("node_layout", ("hex", "rect"))
def test_origin_keyword(node_layout, orientation, n_rows, n_cols):
    """Test setting the origin."""
    graph = TriGraph((n_rows, n_cols))

    assert np.min(graph.x_of_node) == approx(0.0)
    assert np.min(graph.y_of_node) == approx(0.0)

    graph = TriGraph((n_rows, n_cols), origin=(0.5, 0.25))

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
    assert_array_equal(graph.perimeter_nodes, [1, 4, 8, 11, 10, 9, 5, 2, 0])


def test_adjacent_nodes_at_node():
    graph = TriGraph((3, 3), node_layout="hex")
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


# @pytest.mark.skip
@pytest.mark.parametrize("n_cols", (2, 3))
@pytest.mark.parametrize("n_rows", (2, 3))
def test_xy_of_node_rect_vertical(n_rows, n_cols):
    expected = {
        (2, 2): ([0, 1, 0, 1], [0, 0.5, 1, 1.5]),
        (2, 3): ([0, 2, 1, 0, 2, 1], [0, 0, 0.5, 1, 1, 1.5]),
        (3, 2): ([0, 1, 0, 1, 0, 1], [0, 0.5, 1, 1.5, 2, 2.5]),
        (3, 3): ([0, 2, 1, 0, 2, 1, 0, 2, 1], [0, 0, 0.5, 1, 1, 1.5, 2, 2, 2.5]),
    }
    x_of_node = np.empty(n_rows * n_cols, dtype=float)
    y_of_node = np.empty(n_rows * n_cols, dtype=float)
    fill_xy_of_node_rect_vertical((n_rows, n_cols), x_of_node, y_of_node)

    assert np.all(x_of_node == approx(expected[(n_rows, n_cols)][0]))
    assert np.all(y_of_node == approx(expected[(n_rows, n_cols)][1]))


@pytest.mark.parametrize("n_cols", (2, 3))
@pytest.mark.parametrize("n_rows", (1, 2, 3))
def test_xy_of_node_hex_vertical(n_rows, n_cols):
    expected = {
        (1, 2): ([0.5, 0, 0.5], [0, 1, 2]),
        (1, 3): ([0.5, 0, 1, 0.5], [0, 1, 1, 2]),
        (2, 2): ([0.5, 0, 0.5, 0, 0.5], [0, 1, 2, 3, 4]),
        (2, 3): ([0.5, 0, 1, 0.5, 0, 1, 0.5], [0, 1, 1, 2, 3, 3, 4]),
        (3, 2): ([0.5, 0, 0.5, 0, 0.5, 0.0, 0.5], [0, 1, 2, 3, 4, 5, 6]),
        (3, 3): (
            [0.5, 0, 1, 0.5, 0, 1, 0.5, 0, 1, 0.5],
            [0, 1, 1, 2, 3, 3, 4, 5, 5, 6],
        ),
    }
    n_nodes = number_of_nodes(
        (n_rows, n_cols), orientation="vertical", node_layout="hex"
    )
    x_of_node = np.empty(n_nodes, dtype=float)
    y_of_node = np.empty(n_nodes, dtype=float)
    fill_xy_of_node_hex_vertical((n_rows, n_cols), x_of_node, y_of_node)

    assert np.all(x_of_node == approx(expected[(n_rows, n_cols)][0]))
    assert np.all(y_of_node == approx(expected[(n_rows, n_cols)][1]))


@pytest.mark.parametrize("n_cols", (2, 3))
@pytest.mark.parametrize("n_rows", (1, 2, 3))
@pytest.mark.parametrize("orientation", ("horizontal", "vertical"))
@pytest.mark.parametrize("node_layout", ("hex", "rect"))
def test_xy_of_node_lower_left(node_layout, orientation, n_rows, n_cols):
    (x_of_node, y_of_node) = setup_xy_of_node(
        (n_rows, n_cols), orientation=orientation, node_layout=node_layout
    )

    assert np.min(x_of_node) == approx(0.0)
    assert np.min(y_of_node) == approx(0.0)

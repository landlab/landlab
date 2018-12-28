#! /usr/bin/env python
import numpy as np
import pytest

from landlab.graph.hex.perimeternodes import number_of_perimeter_nodes, perimeter_nodes
from landlab.graph.hex.hex import number_of_nodes


@pytest.mark.parametrize("node_layout", ["hex", "rect"])
@pytest.mark.parametrize("orientation", ["vertical", "horizontal"])
@pytest.mark.parametrize("n_cols", np.random.randint(3, high=1000, size=5))
@pytest.mark.parametrize("n_rows", np.random.randint(3, high=1000, size=5))
def test_perimeter_node_count(n_rows, n_cols, orientation, node_layout):
    n_nodes = number_of_perimeter_nodes(
        (n_rows, n_cols), orientation=orientation, node_layout=node_layout
    )
    assert n_nodes > 0
    assert n_nodes < number_of_nodes((n_rows, n_cols), node_layout=node_layout)


EXPECTED_PERIMETER_NODE_COUNT = {
    "rect": {"horizontal": {3: 10, 4: 12}, "vertical": {3: 10, 4: 12}},
    "hex": {"horizontal": {3: 10, 4: 13}, "vertical": {3: 10, 4: 13}},
}


@pytest.mark.parametrize("n_rows", (3, 4))
@pytest.mark.parametrize("orientation", ("horizontal", "vertical"))
@pytest.mark.parametrize("node_layout", ("rect", "hex"))
def test_calc_perimeter_node_count(orientation, node_layout, n_rows):
    if orientation == "horizontal":
        shape = (n_rows, 4)
    else:
        shape = (4, n_rows)
    n_nodes = number_of_perimeter_nodes(
        shape, orientation=orientation, node_layout=node_layout
    )
    assert n_nodes == EXPECTED_PERIMETER_NODE_COUNT[node_layout][orientation][n_rows]


@pytest.mark.parametrize("orientation", ("horizontal", "vertical"))
@pytest.mark.parametrize("n_cols", [3, 7, 12])
@pytest.mark.parametrize("n_rows", [3, 7, 12])
def test_perimeter_nodes_hex(n_rows, n_cols, orientation):
    expected = n_rows * 2 + (n_cols - 2) * 2

    if orientation == "horizontal":
        shape = (n_rows, n_cols)
    else:
        shape = (n_cols, n_rows)

    if shape[0] % 2 == 0:
        expected += 1

    n_nodes = number_of_perimeter_nodes(
        shape, orientation="horizontal", node_layout="hex"
    )
    assert n_nodes == expected


@pytest.mark.parametrize("orientation", ("horizontal", "vertical"))
@pytest.mark.parametrize("n_cols", [3, 7, 12])
@pytest.mark.parametrize("n_rows", [3, 7, 12])
def test_perimeter_nodes_rect(n_rows, n_cols, orientation):
    expected = n_rows * 2 + (n_cols - 2) * 2

    if orientation == "horizontal":
        shape = (n_rows, n_cols)
    else:
        shape = (n_cols, n_rows)
    n_nodes = number_of_perimeter_nodes(
        shape, orientation="horizontal", node_layout="rect"
    )
    assert n_nodes == expected


EXPECTED_PERIMETER_NODES = {
    "rect": {
        "horizontal": {
            3: (3, 7, 11, 10, 9, 8, 4, 0, 1, 2),
            4: (3, 7, 11, 15, 14, 13, 12, 8, 4, 0, 1, 2),
        },
        "vertical": {
            3: (1, 4, 7, 10, 11, 9, 6, 3, 0, 2),
            4: (3, 7, 11, 15, 13, 14, 12, 8, 4, 0, 2, 1),
        },
    },
    "hex": {
        "horizontal": {
            3: (3, 8, 12, 11, 10, 9, 4, 0, 1, 2),
            4: (3, 8, 14, 19, 18, 17, 16, 15, 9, 4, 0, 1, 2),
        },
        "vertical": {
            3: (2, 5, 8, 11, 12, 10, 7, 4, 1, 0),
            4: (2, 6, 10, 14, 18, 19, 17, 15, 11, 7, 3, 1, 0),
        },
    },
}


@pytest.mark.parametrize("n_rows", (3, 4))
@pytest.mark.parametrize("orientation", ("horizontal", "vertical"))
@pytest.mark.parametrize("node_layout", ("rect", "hex"))
def test_calc_perimeter_nodes(orientation, node_layout, n_rows):
    shape = (n_rows, 4)
    if orientation == "vertical":
        shape = (shape[1], shape[0])
    nodes = perimeter_nodes(shape, orientation=orientation, node_layout=node_layout)
    assert tuple(nodes) == EXPECTED_PERIMETER_NODES[node_layout][orientation][n_rows]

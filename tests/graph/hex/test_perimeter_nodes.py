#! /usr/bin/env python
import numpy as np
import pytest
from numpy.testing import assert_array_equal

from landlab.graph.hex.hex import HorizontalHexTriGraph
from landlab.graph.hex.hex import HorizontalRectTriGraph
from landlab.graph.hex.hex import VerticalHexTriGraph
from landlab.graph.hex.hex import VerticalRectTriGraph


@pytest.mark.parametrize(
    "n_cols", (1,) + tuple(np.random.randint(2, high=1000, size=5))
)
@pytest.mark.parametrize(
    "n_rows", (1,) + tuple(np.random.randint(2, high=1000, size=5))
)
def test_perimeter_node_count(n_rows, n_cols, hex_layout):
    n_nodes = hex_layout.number_of_perimeter_nodes((n_rows, n_cols))

    assert n_nodes > 0
    assert n_nodes <= hex_layout.number_of_nodes((n_rows, n_cols))


def test_calc_perimeter_node_count_horizontal_rect():
    assert HorizontalRectTriGraph.number_of_perimeter_nodes((1, 4)) == 4
    assert HorizontalRectTriGraph.number_of_perimeter_nodes((3, 4)) == 10
    assert HorizontalRectTriGraph.number_of_perimeter_nodes((4, 4)) == 12


def test_calc_perimeter_node_count_vertical_rect():
    assert VerticalRectTriGraph.number_of_perimeter_nodes((4, 1)) == 4
    assert VerticalRectTriGraph.number_of_perimeter_nodes((4, 3)) == 10
    assert VerticalRectTriGraph.number_of_perimeter_nodes((4, 4)) == 12


def test_calc_perimeter_node_count_horizontal_hex():
    assert HorizontalHexTriGraph.number_of_perimeter_nodes((1, 4)) == 4
    assert HorizontalHexTriGraph.number_of_perimeter_nodes((3, 4)) == 10
    assert HorizontalHexTriGraph.number_of_perimeter_nodes((4, 4)) == 13


def test_calc_perimeter_node_count_vertical_hex():
    assert VerticalHexTriGraph.number_of_perimeter_nodes((4, 1)) == 4
    assert VerticalHexTriGraph.number_of_perimeter_nodes((4, 3)) == 10
    assert VerticalHexTriGraph.number_of_perimeter_nodes((4, 4)) == 13


@pytest.mark.parametrize("n_cols", [1, 2, 3, 7, 12])
@pytest.mark.parametrize("n_rows", [1, 2, 3, 7, 12])
def test_length_of_perimeter_nodes(hex_layout, n_rows, n_cols):
    assert len(
        hex_layout.perimeter_nodes((n_rows, n_cols))
    ) == hex_layout.number_of_perimeter_nodes((n_rows, n_cols))


def test_calc_perimeter_nodes_horizontal_rect():
    assert_array_equal(
        HorizontalRectTriGraph.perimeter_nodes((3, 4)), (3, 7, 11, 10, 9, 8, 4, 0, 1, 2)
    )
    assert_array_equal(
        HorizontalRectTriGraph.perimeter_nodes((4, 4)),
        (3, 7, 11, 15, 14, 13, 12, 8, 4, 0, 1, 2),
    )

    assert_array_equal(HorizontalRectTriGraph.perimeter_nodes((1, 4)), (3, 2, 1, 0))
    assert_array_equal(HorizontalRectTriGraph.perimeter_nodes((4, 1)), (0, 1, 2, 3))


def test_calc_perimeter_nodes_vertical_rect():
    assert_array_equal(
        VerticalRectTriGraph.perimeter_nodes((4, 3)), (1, 4, 7, 10, 11, 9, 6, 3, 0, 2)
    )
    assert_array_equal(
        VerticalRectTriGraph.perimeter_nodes((4, 4)),
        (3, 7, 11, 15, 13, 14, 12, 8, 4, 0, 2, 1),
    )

    assert_array_equal(VerticalRectTriGraph.perimeter_nodes((4, 1)), (0, 1, 2, 3))
    assert_array_equal(VerticalRectTriGraph.perimeter_nodes((1, 4)), (3, 1, 2, 0))


def test_calc_perimeter_nodes_horizontal_hex():
    assert_array_equal(
        HorizontalHexTriGraph.perimeter_nodes((3, 4)), (8, 12, 11, 10, 9, 4, 0, 1, 2, 3)
    )
    assert_array_equal(
        HorizontalHexTriGraph.perimeter_nodes((4, 4)),
        (14, 19, 18, 17, 16, 15, 9, 4, 0, 1, 2, 3, 8),
    )

    assert_array_equal(HorizontalHexTriGraph.perimeter_nodes((1, 4)), (3, 2, 1, 0))


def test_calc_perimeter_nodes_vertical_hex():
    assert_array_equal(
        VerticalHexTriGraph.perimeter_nodes((4, 3)), (2, 5, 8, 11, 12, 10, 7, 4, 1, 0)
    )
    assert_array_equal(
        VerticalHexTriGraph.perimeter_nodes((4, 4)),
        (2, 6, 10, 14, 18, 19, 17, 15, 11, 7, 3, 1, 0),
    )
    assert_array_equal(
        VerticalHexTriGraph.perimeter_nodes((4, 4)),
        (2, 6, 10, 14, 18, 19, 17, 15, 11, 7, 3, 1, 0),
    )

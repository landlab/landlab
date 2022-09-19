#! /usr/bin/env python
import numpy as np
import pytest
from numpy.testing import assert_array_equal

from landlab.graph.framed_voronoi.framed_voronoi import HorizontalRectVoronoiGraph


@pytest.mark.parametrize(
    "n_cols", (1,) + tuple(np.random.randint(2, high=1000, size=5))
)
@pytest.mark.parametrize(
    "n_rows", (1,) + tuple(np.random.randint(2, high=1000, size=5))
)
def test_rect_perimeter_node_count(n_rows, n_cols, layout_graph):
    n_nodes = layout_graph.number_of_perimeter_nodes((n_rows, n_cols))

    assert n_nodes > 0
    assert n_nodes <= layout_graph.number_of_nodes((n_rows, n_cols))


def test_calc_perimeter_node_count_horizontal_rect():
    assert HorizontalRectVoronoiGraph.number_of_perimeter_nodes((1, 4)) == 4
    assert HorizontalRectVoronoiGraph.number_of_perimeter_nodes((3, 4)) == 10
    assert HorizontalRectVoronoiGraph.number_of_perimeter_nodes((4, 4)) == 12


@pytest.mark.parametrize("n_cols", [1, 2, 3, 7, 12])
@pytest.mark.parametrize("n_rows", [1, 2, 3, 7, 12])
def test_length_of_perimeter_nodes(layout_graph, n_rows, n_cols):
    assert len(
        layout_graph.perimeter_nodes((n_rows, n_cols))
    ) == layout_graph.number_of_perimeter_nodes((n_rows, n_cols))


def test_calc_perimeter_nodes_horizontal_rect():
    assert_array_equal(
        HorizontalRectVoronoiGraph.perimeter_nodes((3, 4)),
        (3, 7, 11, 10, 9, 8, 4, 0, 1, 2),
    )
    assert_array_equal(
        HorizontalRectVoronoiGraph.perimeter_nodes((4, 4)),
        (3, 7, 11, 15, 14, 13, 12, 8, 4, 0, 1, 2),
    )

    assert_array_equal(HorizontalRectVoronoiGraph.perimeter_nodes((1, 4)), (3, 2, 1, 0))
    assert_array_equal(HorizontalRectVoronoiGraph.perimeter_nodes((4, 1)), (0, 1, 2, 3))

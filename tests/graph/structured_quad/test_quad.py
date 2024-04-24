"""Test StructuredQuadGraph."""

import numpy as np
from numpy.testing import assert_array_equal
from pytest import approx
from pytest import mark
from pytest import raises

from landlab.graph import RectilinearGraph
from landlab.graph import StructuredQuadGraph
from landlab.graph import UniformRectilinearGraph
from landlab.graph.structured_quad.structured_quad import StructuredQuadLayoutCython
from landlab.graph.structured_quad.structured_quad import StructuredQuadLayoutPython


def test_graph_is_frozen():
    graph = UniformRectilinearGraph((3, 4))

    assert_array_equal(
        graph.nodes_at_link,
        [
            [0, 1],
            [1, 2],
            [2, 3],
            [0, 4],
            [1, 5],
            [2, 6],
            [3, 7],
            [4, 5],
            [5, 6],
            [6, 7],
            [4, 8],
            [5, 9],
            [6, 10],
            [7, 11],
            [8, 9],
            [9, 10],
            [10, 11],
        ],
    )

    with raises(ValueError):
        graph.nodes_at_link[0] = [1, 0]


def test_graph_can_thaw():
    graph = UniformRectilinearGraph((3, 4))

    assert_array_equal(
        graph.nodes_at_link,
        [
            [0, 1],
            [1, 2],
            [2, 3],
            [0, 4],
            [1, 5],
            [2, 6],
            [3, 7],
            [4, 5],
            [5, 6],
            [6, 7],
            [4, 8],
            [5, 9],
            [6, 10],
            [7, 11],
            [8, 9],
            [9, 10],
            [10, 11],
        ],
    )

    with graph.thawed():
        graph.nodes_at_link[0] = [1, 0]

    assert_array_equal(
        graph.nodes_at_link,
        [
            [1, 0],
            [1, 2],
            [2, 3],
            [0, 4],
            [1, 5],
            [2, 6],
            [3, 7],
            [4, 5],
            [5, 6],
            [6, 7],
            [4, 8],
            [5, 9],
            [6, 10],
            [7, 11],
            [8, 9],
            [9, 10],
            [10, 11],
        ],
    )


@mark.parametrize("layout", (StructuredQuadLayoutCython, StructuredQuadLayoutPython))
def test_layout_links_at_patch(layout):
    assert_array_equal(
        layout.links_at_patch((3, 4)),
        [
            [4, 7, 3, 0],
            [5, 8, 4, 1],
            [6, 9, 5, 2],
            [11, 14, 10, 7],
            [12, 15, 11, 8],
            [13, 16, 12, 9],
        ],
    )


@mark.parametrize(
    "method",
    (
        "links_at_patch",
        "nodes_at_link",
        "horizontal_links",
        "vertical_links",
        "perimeter_nodes",
        "links_at_node",
        "patches_at_link",
        "link_dirs_at_node",
        "patches_at_node",
    ),
)
def test_layouts_match(method):
    assert_array_equal(
        getattr(StructuredQuadLayoutCython, method)((3, 4)),
        getattr(StructuredQuadLayoutPython, method)((3, 4)),
    )


@mark.skip("speed tests")
@mark.parametrize(
    "method",
    (
        "links_at_patch",
        "nodes_at_link",
        "horizontal_links",
        "vertical_links",
        "perimeter_nodes",
        "links_at_node",
        "patches_at_link",
        "link_dirs_at_node",
        "patches_at_node",
    ),
)
@mark.parametrize("size", (10, 11))
def test_layouts_cython_is_faster(method, size):
    from timeit import timeit

    n_rows, n_cols = 3 * 2**size, 4 * 2**size

    def time_method(impl):
        return timeit(
            "{impl}.{method}(({n_rows}, {n_cols}))".format(
                impl=impl, method=method, n_rows=n_rows, n_cols=n_cols
            ),
            setup="from landlab.graph.structured_quad.structured_quad import {impl}".format(
                impl=impl
            ),
            number=1,
        )

    benchmark = time_method("StructuredQuadLayoutPython")
    time = time_method("StructuredQuadLayoutCython")
    speedup = benchmark / time

    assert speedup > 1  # or time < 1e-2


def test_create():
    """Test creating a quad graph."""
    y = [0, 1, 3, 0, 1, 3, 0, 1, 3]
    x = [3, 3, 3, 4, 4, 4, 6, 6, 6]
    graph = StructuredQuadGraph((y, x), shape=(3, 3))

    assert graph.number_of_nodes == 9
    assert graph.number_of_links == 12
    assert graph.number_of_patches == 4


def test_nodes():
    graphs = (
        UniformRectilinearGraph((3, 4)),
        StructuredQuadGraph(
            (
                [[0.0, 0.0, 0.0, 0.0], [1.0, 1.0, 1.0, 1.0], [2.0, 2.0, 2.0, 2.0]],
                [[0.0, 1.0, 2.0, 3.0], [0.0, 1.0, 2.0, 3.0], [0.0, 1.0, 2.0, 3.0]],
            )
        ),
        RectilinearGraph(([0.0, 1.0, 2.0], [0.0, 1.0, 2.0, 3.0])),
    )

    for graph in graphs:
        assert_array_equal(graph.nodes, [[0, 1, 2, 3], [4, 5, 6, 7], [8, 9, 10, 11]])
        with raises(ValueError):
            graph.nodes[0, 0] = 99


def test_perimeter_nodes():
    y = [0, 1, 3, 0, 1, 3, 0, 1, 3]
    x = [3, 3, 3, 4, 4, 4, 6, 6, 6]
    graph = StructuredQuadGraph((y, x), shape=(3, 3))
    assert_array_equal(graph.perimeter_nodes, [2, 5, 8, 7, 6, 3, 0, 1])


def test_length_of_link():
    """Test length of links."""
    y = [0, 1, 3, 0, 1, 3, 0, 1, 3]
    x = [3, 3, 3, 4, 4, 4, 6, 6, 6]
    graph = StructuredQuadGraph((y, x), shape=(3, 3))
    assert graph.length_of_link == approx(
        [1.0, 2.0, 1.0, 1.0, 1.0, 1.0, 2.0, 2.0, 2.0, 2.0, 1.0, 2.0]
    )


def test_area_of_patch():
    """Test areas of patches."""
    y = [0, 0, 0, 1, 1, 1, 3, 3, 3]
    x = [3, 4, 6, 3, 4, 6, 3, 4, 6]
    graph = StructuredQuadGraph((y, x), shape=(3, 3))
    assert graph.area_of_patch == approx([1.0, 2.0, 2.0, 4.0])


def test_nodes_at_patch():
    """Test areas of patches."""
    y = [0, 1, 3, 0, 1, 3, 0, 1, 3]
    x = [3, 3, 3, 4, 4, 4, 6, 6, 6]
    graph = StructuredQuadGraph((y, x), shape=(3, 3), sort=True)
    assert_array_equal(
        graph.nodes_at_patch, [[4, 3, 0, 1], [5, 4, 1, 2], [7, 6, 3, 4], [8, 7, 4, 5]]
    )


def test_patches_at_node():
    """Test areas of patches."""
    y = [0, 1, 3, 0, 1, 3, 0, 1, 3]
    x = [3, 3, 3, 4, 4, 4, 6, 6, 6]
    graph = StructuredQuadGraph((y, x), shape=(3, 3))
    assert_array_equal(
        graph.patches_at_node,
        [
            [0, -1, -1, -1],
            [1, 0, -1, -1],
            [-1, 1, -1, -1],
            [2, -1, -1, 0],
            [3, 2, 0, 1],
            [-1, 3, 1, -1],
            [-1, -1, -1, 2],
            [-1, -1, 2, 3],
            [-1, -1, 3, -1],
        ],
    )


def test_patches_at_link():
    """Test areas of patches."""
    y = [0, 1, 3, 0, 1, 3, 0, 1, 3]
    x = [3, 3, 3, 4, 4, 4, 6, 6, 6]
    graph = StructuredQuadGraph((y, x), shape=(3, 3))
    assert_array_equal(
        graph.patches_at_link,
        [
            [-1, 0],
            [-1, 1],
            [0, -1],
            [1, 0],
            [-1, 1],
            [0, 2],
            [1, 3],
            [2, -1],
            [3, 2],
            [-1, 3],
            [2, -1],
            [3, -1],
        ],
    )


def test_links_at_patch():
    """Test areas of patches."""
    y = [0, 1, 3, 0, 1, 3, 0, 1, 3]
    x = [3, 3, 3, 4, 4, 4, 6, 6, 6]
    graph = StructuredQuadGraph((y, x), shape=(3, 3))
    assert_array_equal(
        graph.links_at_patch, [[3, 5, 2, 0], [4, 6, 3, 1], [8, 10, 7, 5], [9, 11, 8, 6]]
    )


def test_nodes_at_link():
    """Test areas of patches."""
    y = [0, 1, 3, 0, 1, 3, 0, 1, 3]
    x = [3, 3, 3, 4, 4, 4, 6, 6, 6]
    graph = StructuredQuadGraph((y, x), shape=(3, 3))
    assert_array_equal(
        graph.nodes_at_link,
        [
            [0, 1],
            [1, 2],
            [0, 3],
            [1, 4],
            [2, 5],
            [3, 4],
            [4, 5],
            [3, 6],
            [4, 7],
            [5, 8],
            [6, 7],
            [7, 8],
        ],
    )


def test_links_at_node():
    """Test areas of patches."""
    y = [0, 1, 3, 0, 1, 3, 0, 1, 3]
    x = [3, 3, 3, 4, 4, 4, 6, 6, 6]
    graph = StructuredQuadGraph((y, x), shape=(3, 3))
    assert_array_equal(
        graph.links_at_node,
        [
            [0, 2, -1, -1],
            [1, 3, 0, -1],
            [-1, 4, 1, -1],
            [5, 7, -1, 2],
            [6, 8, 5, 3],
            [-1, 9, 6, 4],
            [10, -1, -1, 7],
            [11, -1, 10, 8],
            [-1, -1, 11, 9],
        ],
    )


def test_link_dirs_at_node():
    """Test areas of patches."""
    y = [0, 1, 3, 0, 1, 3, 0, 1, 3]
    x = [3, 3, 3, 4, 4, 4, 6, 6, 6]
    graph = StructuredQuadGraph((y, x), shape=(3, 3))
    assert_array_equal(
        graph.link_dirs_at_node,
        [
            [-1, -1, 0, 0],
            [-1, -1, 1, 0],
            [0, -1, 1, 0],
            [-1, -1, 0, 1],
            [-1, -1, 1, 1],
            [0, -1, 1, 1],
            [-1, 0, 0, 1],
            [-1, 0, 1, 1],
            [0, 0, 1, 1],
        ],
    )


def test_link_dirs_at_node_raster():
    graph = UniformRectilinearGraph((4, 3))
    assert_array_equal(
        graph.link_dirs_at_node,
        [
            [-1, -1, 0, 0],
            [-1, -1, 1, 0],
            [0, -1, 1, 0],
            [-1, -1, 0, 1],
            [-1, -1, 1, 1],
            [0, -1, 1, 1],
            [-1, -1, 0, 1],
            [-1, -1, 1, 1],
            [0, -1, 1, 1],
            [-1, 0, 0, 1],
            [-1, 0, 1, 1],
            [0, 0, 1, 1],
        ],
    )
    assert graph.link_dirs_at_node.dtype == np.int8

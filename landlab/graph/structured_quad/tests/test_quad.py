"""Test StructuredQuadGraph."""
from numpy.testing import assert_array_almost_equal, assert_array_equal

from landlab.graph import StructuredQuadGraph


def test_create():
    """Test creating a quad graph."""
    y = [0, 1, 3, 0, 1, 3, 0, 1, 3]
    x = [3, 3, 3, 4, 4, 4, 6, 6, 6]
    graph = StructuredQuadGraph((y, x), shape=(3, 3))

    assert graph.number_of_nodes == 9
    assert graph.number_of_links == 12
    assert graph.number_of_patches == 4


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
    assert_array_almost_equal(
        graph.length_of_link,
        [1.0, 2.0, 1.0, 1.0, 1.0, 1.0, 2.0, 2.0, 2.0, 2.0, 1.0, 2.0],
    )


def test_area_of_patch():
    """Test areas of patches."""
    y = [0, 1, 3, 0, 1, 3, 0, 1, 3]
    x = [3, 3, 3, 4, 4, 4, 6, 6, 6]
    graph = StructuredQuadGraph((y, x), shape=(3, 3))
    assert_array_almost_equal(graph.area_of_patch, [1.0, 2.0, 2.0, 4.0])


def test_nodes_at_patch():
    """Test areas of patches."""
    y = [0, 1, 3, 0, 1, 3, 0, 1, 3]
    x = [3, 3, 3, 4, 4, 4, 6, 6, 6]
    graph = StructuredQuadGraph((y, x), shape=(3, 3))
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

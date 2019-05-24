from numpy.testing import assert_array_almost_equal, assert_array_equal

from landlab.graph import Graph

r"""
For these tests the nodes are given column-by-column::

    (1) -- (3) --- (4)
     |      |     /
    (0)     |   /
        \   | /
           (2)

Once sorted, the node numbering becomes::

    (2) -- (3) --- (4)
     |      |     /
    (1)     |   /
        \   | /
           (0)
"""
NODE_X = (0, 0, 1, 1, 2, 2)
NODE_Y = (0, 1, 0, 1, 1, 0)
NODES_AT_LINK = ((0, 1), (0, 2), (1, 3), (2, 3), (2, 5), (3, 4), (4, 5))
LINKS_AT_PATCH = ((3, 2, 0, 1), (3, 5, 6, 4))

# NODE_X = (0, 0, 1, 1, 2)
# NODE_Y = (1, 2, 0, 2, 2)
# NODES_AT_LINK = ((0, 1), (0, 2), (1, 3), (2, 3), (2, 4), (3, 4))
# LINKS_AT_PATCH = ((3, 2, 0, 1), (3, 5, 4))
# LINKS_AT_PATCH = ([3, 2, 0, 1, 3, 5, 4], [4, 3])


def test_create_graph_with_nodes():
    """Create a graph of unconnected nodes."""
    graph = Graph((NODE_Y, NODE_X))

    assert_array_almost_equal(graph.x_of_node, [0.0, 1.0, 2.0, 0.0, 1.0, 2.0])
    assert_array_almost_equal(graph.y_of_node, [0.0, 0.0, 0.0, 1.0, 1.0, 1.0])


def test_create_graph_with_links():
    """Create a graph of connected nodes."""
    graph = Graph((NODE_Y, NODE_X), links=NODES_AT_LINK)

    assert_array_almost_equal(graph.x_of_node, [0.0, 1.0, 2.0, 0.0, 1.0, 2.0])
    assert_array_almost_equal(graph.y_of_node, [0.0, 0.0, 0.0, 1.0, 1.0, 1.0])
    assert_array_almost_equal(
        graph.nodes_at_link, [[0, 1], [1, 2], [0, 3], [1, 4], [2, 5], [3, 4], [4, 5]]
    )


def test_graph_link_heads():
    """Test nodes at link heads."""
    graph = Graph((NODE_Y, NODE_X), links=NODES_AT_LINK)

    assert_array_equal(graph.node_at_link_head, [1, 2, 3, 4, 5, 4, 5])


def test_graph_link_tails():
    """Test nodes at link tails."""
    graph = Graph((NODE_Y, NODE_X), links=NODES_AT_LINK)

    assert_array_equal(graph.node_at_link_tail, [0, 1, 0, 1, 2, 3, 4])


def test_graph_links_at_node():
    """Test links at nodes without rotational sorting."""
    graph = Graph((NODE_Y, NODE_X), links=NODES_AT_LINK)

    assert_array_equal(
        graph.links_at_node,
        [[0, 2, -1], [1, 3, 0], [4, 1, -1], [5, 2, -1], [6, 5, 3], [6, 4, -1]],
    )


def test_graph_link_dirs_at_node():
    """Test links directions at nodes without rotational sorting."""
    graph = Graph((NODE_Y, NODE_X), links=NODES_AT_LINK)

    assert_array_equal(
        graph.link_dirs_at_node,
        [[-1, -1, 0], [-1, -1, 1], [-1, 1, 0], [-1, 1, 0], [-1, 1, 1], [1, 1, 0]],
    )


def test_links_at_patch_ccw():
    """Test links at patch with rotational sorting."""
    graph = Graph((NODE_Y, NODE_X), links=NODES_AT_LINK, patches=LINKS_AT_PATCH)

    assert_array_equal(graph.links_at_patch, [[3, 5, 2, 0], [4, 6, 3, 1]])


def test_nodes_at_patch_ccw():
    """Test nodes at patch with rotational sorting."""
    graph = Graph((NODE_Y, NODE_X), links=NODES_AT_LINK, patches=LINKS_AT_PATCH)

    assert_array_equal(graph.nodes_at_patch, [[4, 3, 0, 1], [5, 4, 1, 2]])

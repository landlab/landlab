from numpy.testing import assert_array_equal, assert_array_almost_equal

from landlab.graph import Graph


"""
The nodes are given column-by-column::

    (1) -- (3) -- (5)
     |      |      |
    (0) -- (2) -- (4)

"""

NODE_X = (0, 0, 1, 1, 2, 2)
NODE_Y = (0, 1, 0, 1, 0, 1)
NODES_AT_LINK = ((0, 1), (0, 2), (1, 3),
                 (2, 3), (2, 4), (3, 5),
                 (4, 5))
LINKS_AT_PATCH = ((3, 2, 0, 1), (6, 5, 3, 4))


def test_create_graph_with_nodes():
    """Create a graph of unconnected nodes."""
    graph = Graph((NODE_Y, NODE_X))

    assert_array_almost_equal(graph.x_of_node, NODE_X)
    assert_array_almost_equal(graph.y_of_node, NODE_Y)


def test_create_graph_with_nodes_xy_sort():
    """Create a graph with nodes sort by x, then y."""
    graph = Graph((NODE_Y, NODE_X), xy_sort=True)

    assert_array_almost_equal(graph.x_of_node, [0., 1., 2., 0., 1., 2.])
    assert_array_almost_equal(graph.y_of_node, [0., 0., 0., 1., 1., 1.])


def test_create_graph_with_links():
    """Create a graph of connected nodes."""
    graph = Graph((NODE_Y, NODE_X), links=NODES_AT_LINK)

    assert_array_almost_equal(graph.x_of_node, NODE_X)
    assert_array_almost_equal(graph.y_of_node, NODE_Y)
    assert_array_almost_equal(graph.nodes_at_link, NODES_AT_LINK)


def test_create_graph_with_links_xy_sort():
    """Create a graph of connected nodes, sort links by x, then y."""
    graph = Graph((NODE_Y, NODE_X), links=NODES_AT_LINK, xy_sort=True)

    assert_array_almost_equal(graph.x_of_node, [0., 1., 2., 0., 1., 2.])
    assert_array_almost_equal(graph.y_of_node, [0., 0., 0., 1., 1., 1.])
    assert_array_almost_equal(graph.nodes_at_link,
                             [(0, 1), (1, 2),
                              (0, 3), (1, 4), (2, 5),
                              (3, 4), (4, 5)])


def test_create_graph_with_links_ccw():
    """Create a graph of connected nodes."""
    graph = Graph((NODE_Y, NODE_X), links=NODES_AT_LINK, rot_sort=True)

    assert_array_almost_equal(graph.x_of_node, NODE_X)
    assert_array_almost_equal(graph.y_of_node, NODE_Y)
    assert_array_almost_equal(graph.nodes_at_link, NODES_AT_LINK)
    assert_array_almost_equal(graph.links_at_node,
                              [[1, 0, -1], [2, 0, -1],
                               [4, 3,  1], [5, 2,  3],
                               [6, 4, -1], [5, 6, -1]])


def test_graph_link_heads():
    graph = Graph((NODE_Y, NODE_X), links=NODES_AT_LINK, rot_sort=True)

    assert_array_equal(graph.node_at_link_head, [1, 2, 3, 3, 4, 5, 5])


def test_graph_link_tails():
    graph = Graph((NODE_Y, NODE_X), links=NODES_AT_LINK, rot_sort=True)

    assert_array_equal(graph.node_at_link_tail,[0, 0, 1, 2, 2, 3, 4])


def test_graph_links_at_node():
    graph = Graph((NODE_Y, NODE_X), links=NODES_AT_LINK, rot_sort=True)

    assert_array_almost_equal(graph.links_at_node,
                              [[1, 0, -1], [2, 0, -1],
                               [4, 3,  1], [5, 2,  3],
                               [6, 4, -1], [5, 6, -1]])


def test_graph_link_dirs_at_node():
    graph = Graph((NODE_Y, NODE_X), links=NODES_AT_LINK, rot_sort=True)

    assert_array_equal(graph.link_dirs_at_node,
                       [[-1, -1, 0], [-1, 1, 0],
                        [-1, -1, 1], [-1, 1, 1],
                        [-1,  1, 0], [ 1, 1, 0]])


def test_create_graph_with_patches():
    graph = Graph((NODE_Y, NODE_X), links=NODES_AT_LINK,
                  patches=LINKS_AT_PATCH)

    assert_array_equal(graph.links_at_patch, LINKS_AT_PATCH)


def test_graph_nodes_at_patch():
    graph = Graph((NODE_Y, NODE_X), links=NODES_AT_LINK,
                  patches=LINKS_AT_PATCH)

    assert_array_equal(graph.nodes_at_patch,
                       [[0, 1, 2, 3], [2, 3, 4, 5]])

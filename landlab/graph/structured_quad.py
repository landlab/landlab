import numpy as np

from .unstructured import Graph


def setup_links_at_patch(shape):
    """Get links that define patches for a raster grid.

    Examples
    --------
    >>> from landlab.graph.structured_quad import setup_links_at_patch
    >>> setup_links_at_patch((3, 4)) # doctest: +NORMALIZE_WHITESPACE
    array([[ 4,  7,  3,  0], [ 5,  8,  4,  1], [ 6,  9,  5,  2],
           [11, 14, 10,  7], [12, 15, 11,  8], [13, 16, 12,  9]])
    """
    links = np.arange((2 * shape[1] - 1) * shape[0])
    links[-shape[1]:] = -1
    links.shape = (shape[0], 2 * shape[1] - 1)
    horizontal_links = links[:, :shape[1] - 1].reshape(shape[0], shape[1] - 1)
    vertical_links = links[:-1, - shape[1]:].reshape(shape[0] - 1, shape[1])
    links_at_patch = np.vstack((vertical_links[:, 1:].flat,
                                horizontal_links[1:, :].flat,
                                vertical_links[:, :-1].flat,
                                horizontal_links[:-1, :].flat)).T

    return links_at_patch


def setup_nodes_at_link(shape):
    """
    Examples
    --------
    >>> from landlab.graph.structured_quad import setup_nodes_at_link
    >>> setup_nodes_at_link((3, 4)) # doctest: +NORMALIZE_WHITESPACE
    array([[ 0,  1], [ 1,  2], [ 2,  3],
           [ 0,  4], [ 1,  5], [ 2,  6], [ 3,  7],
           [ 4,  5], [ 5,  6], [ 6,  7],
           [ 4,  8], [ 5,  9], [ 6, 10], [ 7, 11],
           [ 8,  9], [ 9, 10], [10, 11]])
    """
    nodes = np.arange((shape[0] + 1) * shape[1]).reshape((shape[0] + 1,
                                                          shape[1]))
    nodes[-1, :] = -1

    link_tails = np.hstack((nodes[:-1, :-1], nodes[:-1, :]))
    link_heads = np.hstack((nodes[:-1, 1:], nodes[1:, :]))

    nodes_at_link = np.vstack((link_tails.flat, link_heads.flat))

    return nodes_at_link.T[:- (shape[1])]


class RasterGraph(Graph):

    """Graph of a structured grid of quadrilaterals.

    Examples
    --------
    >>> from landlab.graph import RasterGraph
    >>> graph = RasterGraph((4, 3))
    >>> graph.number_of_nodes
    12
    >>> graph.links_at_node
    array([[ 0,  2, -1, -1], [ 0,  1,  3, -1], [ 1,  4, -1, -1],
           [ 2,  5,  7, -1], [ 3,  5,  6,  8], [ 4,  6,  9, -1],
           [ 7, 10, 12, -1], [ 8, 10, 11, 13], [ 9, 11, 14, -1],
           [12, 15, -1, -1], [13, 15, 16, -1], [14, 16, -1, -1]])
    >>> graph.nodes_at_link # doctest: +NORMALIZE_WHITESPACE
    array([[ 0,  1], [ 1,  2], [ 0,  3], [ 1,  4], [ 2,  5],
           [ 3,  4], [ 4,  5], [ 3,  6], [ 4,  7], [ 5,  8],
           [ 6,  7], [ 7,  8], [ 6,  9], [ 7, 10], [ 8, 11],
           [ 9, 10], [10, 11]])
    >>> graph.links_at_patch # doctest: +NORMALIZE_WHITESPACE
    array([[ 3,  5,  2,  0], [ 4,  6,  3,  1],
           [ 8, 10,  7,  5], [ 9, 11,  8,  6],
           [13, 15, 12, 10], [14, 16, 13, 11]])
    >>> graph.nodes_at_patch # doctest: +NORMALIZE_WHITESPACE
    array([[ 0,  1,  3,  4], [ 1,  2,  4,  5],
           [ 3,  4,  6,  7], [ 4,  5,  7,  8],
           [ 6,  7,  9, 10], [ 7,  8, 10, 11]])
    """

    def __init__(self, shape, spacing=(1., 1.), origin=(0., 0.)):

        rows = np.arange(shape[0]) * spacing[0] + origin[0]
        cols = np.arange(shape[1]) * spacing[1] + origin[1]
        node_y, node_x = np.meshgrid(rows, cols)
        # nodes = np.arange(shape[0] * shape[1]).reshape(shape)

        nodes_at_link = setup_nodes_at_link(shape)
        links_at_patch = setup_links_at_patch(shape)

        super(RasterGraph, self).__init__((node_y.flat, node_x.flat),
                                          links=nodes_at_link,
                                          patches=links_at_patch)

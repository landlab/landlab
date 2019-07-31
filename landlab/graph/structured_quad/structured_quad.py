from functools import lru_cache
import numpy as np

from ...utils.decorators import read_only_array, store_result_in_grid
from ..graph import Graph
from .ext.at_link import fill_nodes_at_link, fill_patches_at_link
from .ext.at_node import (
    fill_link_dirs_at_node,
    fill_links_at_node,
    fill_patches_at_node,
    fill_perimeter_nodes,
)
from .ext.at_patch import fill_links_at_patch


def setup_links_at_patch(shape):
    """Get links that define patches for a raster grid.

    Examples
    --------
    >>> from landlab.graph.structured_quad.structured_quad import setup_links_at_patch
    >>> setup_links_at_patch((3, 4)) # doctest: +NORMALIZE_WHITESPACE
    array([[ 4,  7,  3,  0], [ 5,  8,  4,  1], [ 6,  9,  5,  2],
           [11, 14, 10,  7], [12, 15, 11,  8], [13, 16, 12,  9]])
    """
    n_patches = (shape[0] - 1) * (shape[1] - 1)
    links_at_patch = np.empty((n_patches, 4), dtype=int)
    fill_links_at_patch(shape, links_at_patch)
    return links_at_patch


def setup_nodes_at_link(shape):
    """
    Examples
    --------
    >>> from landlab.graph.structured_quad.structured_quad import setup_nodes_at_link
    >>> setup_nodes_at_link((3, 4)) # doctest: +NORMALIZE_WHITESPACE
    array([[ 0,  1], [ 1,  2], [ 2,  3],
           [ 0,  4], [ 1,  5], [ 2,  6], [ 3,  7],
           [ 4,  5], [ 5,  6], [ 6,  7],
           [ 4,  8], [ 5,  9], [ 6, 10], [ 7, 11],
           [ 8,  9], [ 9, 10], [10, 11]])
    """
    n_links = shape[0] * (shape[1] - 1) + (shape[0] - 1) * shape[1]
    nodes_at_link = np.empty((n_links, 2), dtype=int)
    fill_nodes_at_link(shape, nodes_at_link)

    return nodes_at_link


class StructuredQuadGraphTopology:
    def __init__(self, shape):
        self._shape = tuple(shape)

    @property
    def shape(self):
        return self._shape

    @property
    def number_of_node_rows(self):
        return self._shape[0]

    @property
    def number_of_node_columns(self):
        return self._shape[1]

    @property
    @lru_cache()
    def nodes_at_right_edge(self):
        return np.arange(self.shape[1] - 1, np.prod(self.shape), self.shape[1])

    @property
    @lru_cache()
    def nodes_at_top_edge(self):
        return np.arange(self.number_of_nodes - self.shape[1], np.prod(self.shape))

    @property
    @lru_cache()
    def nodes_at_left_edge(self):
        return np.arange(0, np.prod(self.shape), self.shape[1])

    @property
    @lru_cache()
    def nodes_at_bottom_edge(self):
        return np.arange(self.shape[1])

    def nodes_at_edge(self, edge):
        if edge not in ("right", "top", "left", "bottom"):
            raise ValueError("value for edge not understood")
        return getattr(self, "nodes_at_{edge}_edge".format(edge=edge))

    @property
    def nodes_at_corners_of_grid(self):
        """Nodes at corners of grid.

        The nodes at at the corners of the grid. The nodes are returned
        counterclockwise starting with the upper-right.

        Return
        ------
        tuple of int
            Nodes at the four corners.

        Examples
        --------
        >>> from landlab.graph import UniformRectilinearGraph
        >>> graph = UniformRectilinearGraph((4, 5))
        >>> graph.nodes_at_corners_of_grid
        (19, 15, 0, 4)
        """
        return (
            self.number_of_nodes - 1,
            self.number_of_nodes - self.number_of_node_columns,
            0,
            self.number_of_node_columns - 1,
        )

    @lru_cache()
    def horizontal_links(self):
        from .ext.at_link import fill_horizontal_links

        n_horizontal_links = self.shape[0] * (self.shape[1] - 1)
        horizontal_links = np.empty(n_horizontal_links, dtype=int)
        fill_horizontal_links(self.shape, horizontal_links)

        return horizontal_links

    @property
    @lru_cache()
    def vertical_links(self):
        from .ext.at_link import fill_vertical_links

        n_vertical_links = (self.shape[0] - 1) * self.shape[1]
        vertical_links = np.empty(n_vertical_links, dtype=int)
        fill_vertical_links(self.shape, vertical_links)

        return vertical_links

    @property
    @lru_cache()
    def perimeter_nodes(self):
        n_perimeter_nodes = 2 * self.shape[0] + 2 * (self.shape[1] - 2)
        perimeter_nodes = np.empty(n_perimeter_nodes, dtype=int)
        fill_perimeter_nodes(self.shape, perimeter_nodes)

        return perimeter_nodes

    @property
    @lru_cache()
    def links_at_node(self):
        n_nodes = self.shape[0] * self.shape[1]
        links_at_node = np.empty((n_nodes, 4), dtype=int)
        fill_links_at_node(self.shape, links_at_node)

        return links_at_node

    @property
    @lru_cache()
    def link_dirs_at_node(self):
        n_nodes = self.shape[0] * self.shape[1]
        link_dirs_at_node = np.empty((n_nodes, 4), dtype=np.int8)
        fill_link_dirs_at_node(self.shape, link_dirs_at_node)
        return link_dirs_at_node

    @property
    @lru_cache()
    def patches_at_node(self):
        n_nodes = self.shape[0] * self.shape[1]
        patches_at_node = np.empty((n_nodes, 4), dtype=int)
        fill_patches_at_node(self.shape, patches_at_node)

        return patches_at_node

    @property
    @lru_cache()
    @read_only_array
    def patches_at_link(self):
        n_links = self.shape[0] * (self.shape[1] - 1) + (self.shape[0] - 1) * self.shape[1]
        patches_at_link = np.empty((n_links, 2), dtype=int)
        fill_patches_at_link(self.shape, patches_at_link)

        return patches_at_link


class StructuredQuadGraphExtras(StructuredQuadGraphTopology, Graph):

    def __init__(self, node_y_and_x):
        StructuredQuadGraphTopology.__init__(self, node_y_and_x[0].shape)
        Graph.__init__(
            self,
            node_y_and_x,
            links=setup_nodes_at_link(self.shape),
            patches=setup_links_at_patch(self.shape),
            sort=False,
        )


class StructuredQuadGraph(StructuredQuadGraphExtras):
    def __init__(self, coords, shape=None):
        node_y, node_x = (
            np.asarray(coords[0], dtype=float),
            np.asarray(coords[1], dtype=float),
        )

        if shape:
            node_y.shape = shape
            node_x.shape = shape
        else:
            node_x.shape = node_y.shape

        if node_y.shape != node_x.shape:
            raise ValueError("shape mismatch in node x and y coordinates")

        StructuredQuadGraphExtras.__init__(self, (node_y, node_x))

    @staticmethod
    def setup_node_y_and_x(yx_at_node, shape=None):
        node_y = np.asarray(yx_at_node[0], dtype=float)
        node_x = np.asarray(yx_at_node[1], dtype=float)

        if shape:
            node_y.shape = shape
            node_x.shape = shape
        else:
            node_x.shape = node_y.shape

        if node_y.shape != node_x.shape:
            raise ValueError("shape mismatch in node x and y coordinates")

        return (node_y, node_x)


class RectilinearGraph(StructuredQuadGraphExtras):

    """Graph of a rectlinear grid of nodes.

    Examples
    --------
    >>> from landlab.graph import RectilinearGraph
    >>> graph = RectilinearGraph(([0, 1, 2, 3], [1, 4, 8]))
    >>> graph.number_of_nodes
    12
    >>> graph.y_of_node # doctest: +NORMALIZE_WHITESPACE
    array([ 0., 0., 0.,
            1., 1., 1.,
            2., 2., 2.,
            3., 3., 3.])
    >>> graph.x_of_node # doctest: +NORMALIZE_WHITESPACE
    array([ 1., 4., 8.,
            1., 4., 8.,
            1., 4., 8.,
            1., 4., 8.])
    """
    def __init__(self, nodes):
        rows = np.asarray(nodes[0], dtype=float)
        cols = np.asarray(nodes[1], dtype=float)
        node_y_and_x = np.meshgrid(rows, cols, indexing="ij")

        StructuredQuadGraphExtras.__init__(self, node_y_and_x)

    @staticmethod
    def setup_node_y_and_x(coords):
        rows = np.asarray(coords[0], dtype=float)
        cols = np.asarray(coords[1], dtype=float)

        return np.meshgrid(rows, cols, indexing="ij")


class UniformRectilinearGraph(StructuredQuadGraphExtras):

    """Graph of a structured grid of quadrilaterals.

    Examples
    --------
    >>> from landlab.graph import UniformRectilinearGraph
    >>> graph = UniformRectilinearGraph((4, 3), spacing=(1, 2), origin=(-1, 0))
    >>> graph.number_of_nodes
    12
    >>> graph.y_of_node # doctest: +NORMALIZE_WHITESPACE
    array([-1., -1., -1.,
            0.,  0.,  0.,
            1.,  1.,  1.,
            2.,  2.,  2.])
    >>> graph.x_of_node # doctest: +NORMALIZE_WHITESPACE
    array([ 0.,  2.,  4.,
            0.,  2.,  4.,
            0.,  2.,  4.,
            0.,  2.,  4.])
    >>> graph.links_at_node # doctest: +NORMALIZE_WHITESPACE
    array([[ 0,  2, -1, -1], [ 1,  3,  0, -1], [-1,  4,  1, -1],
           [ 5,  7, -1,  2], [ 6,  8,  5,  3], [-1,  9,  6,  4],
           [10, 12, -1,  7], [11, 13, 10,  8], [-1, 14, 11,  9],
           [15, -1, -1, 12], [16, -1, 15, 13], [-1, -1, 16, 14]])
    >>> graph.link_dirs_at_node # doctest: +NORMALIZE_WHITESPACE
    array([[-1, -1,  0,  0], [-1, -1,  1,  0], [ 0, -1,  1,  0],
           [-1, -1,  0,  1], [-1, -1,  1,  1], [ 0, -1,  1,  1],
           [-1, -1,  0,  1], [-1, -1,  1,  1], [ 0, -1,  1,  1],
           [-1,  0,  0,  1], [-1,  0,  1,  1], [ 0,  0,  1,  1]], dtype=int8)
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
    array([[ 4,  3,  0,  1], [ 5,  4,  1,  2],
           [ 7,  6,  3,  4], [ 8,  7,  4,  5],
           [10,  9,  6,  7], [11, 10,  7,  8]])
    """

    def __init__(self, shape, spacing=1.0, origin=0.0):
        spacing = np.broadcast_to(spacing, 2)
        origin = np.broadcast_to(origin, 2)

        rows = np.arange(shape[0], dtype=float) * spacing[0] + origin[0]
        cols = np.arange(shape[1], dtype=float) * spacing[1] + origin[1]

        node_y_and_x = np.meshgrid(rows, cols, indexing="ij")

        StructuredQuadGraphExtras.__init__(self, node_y_and_x)

        self._spacing = tuple(spacing)
        self._origin = tuple(origin)

    @property
    def spacing(self):
        return self._spacing

    @property
    def origin(self):
        return self._origin

    @property
    def dx(self):
        return self._spacing[1]

    @property
    def dy(self):
        return self._spacing[0]

from abc import ABC, abstractmethod
from functools import lru_cache

import numpy as np

from ...utils.decorators import read_only_array
from ..graph import Graph


class StructuredQuadLayout(ABC):
    @staticmethod
    def corner_nodes(shape):
        n_rows, n_cols = shape
        return (n_rows * n_cols - 1, (n_rows - 1) * n_cols, 0, n_cols - 1)

    @abstractmethod
    def links_at_patch(shape):
        ...

    @abstractmethod
    def nodes_at_link(shape):
        ...

    @abstractmethod
    def horizontal_links(shape):
        ...

    @abstractmethod
    def vertical_links(shape):
        ...

    @abstractmethod
    def perimeter_nodes(shape):
        ...

    @abstractmethod
    def links_at_node(shape):
        ...

    @abstractmethod
    def patches_at_link(shape):
        ...

    @abstractmethod
    def link_dirs_at_node(shape):
        ...

    @abstractmethod
    def patches_at_node(shape):
        ...


class StructuredQuadLayoutCython(StructuredQuadLayout):
    @staticmethod
    def links_at_patch(shape):
        """Get links that define patches for a raster grid.

        Examples
        --------
        >>> from landlab.graph.structured_quad.structured_quad import StructuredQuadLayoutCython
        >>> StructuredQuadLayoutCython.links_at_patch((3, 4))
        array([[ 4,  7,  3,  0], [ 5,  8,  4,  1], [ 6,  9,  5,  2],
               [11, 14, 10,  7], [12, 15, 11,  8], [13, 16, 12,  9]])
        """
        from .ext.at_patch import fill_links_at_patch

        n_patches = (shape[0] - 1) * (shape[1] - 1)
        links_at_patch = np.empty((n_patches, 4), dtype=int)
        fill_links_at_patch(shape, links_at_patch)
        return links_at_patch

    @staticmethod
    def nodes_at_link(shape):
        """
        Examples
        --------
        >>> from landlab.graph.structured_quad.structured_quad import StructuredQuadLayoutCython
        >>> StructuredQuadLayoutCython.nodes_at_link((3, 4))
        array([[ 0,  1], [ 1,  2], [ 2,  3],
               [ 0,  4], [ 1,  5], [ 2,  6], [ 3,  7],
               [ 4,  5], [ 5,  6], [ 6,  7],
               [ 4,  8], [ 5,  9], [ 6, 10], [ 7, 11],
               [ 8,  9], [ 9, 10], [10, 11]])
        """
        from .ext.at_link import fill_nodes_at_link

        n_links = shape[0] * (shape[1] - 1) + (shape[0] - 1) * shape[1]
        nodes_at_link = np.empty((n_links, 2), dtype=int)
        fill_nodes_at_link(shape, nodes_at_link)

        return nodes_at_link

    @staticmethod
    def horizontal_links(shape):
        from .ext.at_link import fill_horizontal_links

        n_horizontal_links = shape[0] * (shape[1] - 1)
        horizontal_links = np.empty(n_horizontal_links, dtype=int)
        fill_horizontal_links(shape, horizontal_links)

        return horizontal_links

    @staticmethod
    def vertical_links(shape):
        from .ext.at_link import fill_vertical_links

        n_vertical_links = (shape[0] - 1) * shape[1]
        vertical_links = np.empty(n_vertical_links, dtype=int)
        fill_vertical_links(shape, vertical_links)

        return vertical_links

    @staticmethod
    def perimeter_nodes(shape):
        from .ext.at_node import fill_perimeter_nodes

        n_perimeter_nodes = 2 * shape[0] + 2 * (shape[1] - 2)
        perimeter_nodes = np.empty(n_perimeter_nodes, dtype=int)
        fill_perimeter_nodes(shape, perimeter_nodes)

        return perimeter_nodes

    @staticmethod
    def links_at_node(shape):
        from .ext.at_node import fill_links_at_node

        n_nodes = shape[0] * shape[1]
        links_at_node = np.empty((n_nodes, 4), dtype=int)
        fill_links_at_node(shape, links_at_node)

        return links_at_node

    @staticmethod
    def patches_at_link(shape):
        from .ext.at_link import fill_patches_at_link

        n_links = shape[0] * (shape[1] - 1) + (shape[0] - 1) * shape[1]
        patches_at_link = np.empty((n_links, 2), dtype=int)
        fill_patches_at_link(shape, patches_at_link)

        return patches_at_link

    @staticmethod
    def link_dirs_at_node(shape):
        from .ext.at_node import fill_link_dirs_at_node

        n_nodes = shape[0] * shape[1]
        link_dirs_at_node = np.empty((n_nodes, 4), dtype=np.int8)
        fill_link_dirs_at_node(shape, link_dirs_at_node)
        return link_dirs_at_node

    @staticmethod
    def patches_at_node(shape):
        from .ext.at_node import fill_patches_at_node

        n_nodes = shape[0] * shape[1]
        patches_at_node = np.empty((n_nodes, 4), dtype=int)
        fill_patches_at_node(shape, patches_at_node)

        return patches_at_node


class StructuredQuadLayoutPython(StructuredQuadLayout):
    @staticmethod
    def links_at_patch(shape):
        n_rows, n_cols = shape
        n_patches = (shape[0] - 1) * (shape[1] - 1)
        links_at_patch = np.empty((4, n_patches), dtype=int)

        patches = np.arange(n_patches, dtype=int).reshape((n_rows - 1, n_cols - 1))
        south_links = patches + np.arange(n_rows - 1).reshape((n_rows - 1, 1)) * n_cols
        links_at_patch[3, :] = south_links.flat
        links_at_patch[2, :] = links_at_patch[3, :] + n_cols - 1
        links_at_patch[1, :] = links_at_patch[3, :] + 2 * n_cols - 1
        links_at_patch[0, :] = links_at_patch[2, :] + 1

        return links_at_patch.T

    @staticmethod
    def nodes_at_link(shape):
        n_rows, n_cols = shape

        nodes_at_link = np.empty(
            (2, n_rows * (n_cols - 1) + (n_rows - 1) * n_cols), dtype=int
        )
        nodes = np.arange(n_rows * n_cols, dtype=int).reshape((n_rows, n_cols))

        nodes_at_link[0, -(n_cols - 1) :] = nodes[-1, :-1]
        nodes_at_link[1, -(n_cols - 1) :] = nodes[-1, 1:]

        head_nodes = nodes_at_link[0, : -(n_cols - 1)].reshape(
            (n_rows - 1, 2 * n_cols - 1)
        )
        head_nodes[:, : n_cols - 1] = nodes[:-1, :-1]
        head_nodes[:, n_cols - 1 :] = nodes[:-1, :]

        tail_nodes = nodes_at_link[1, : -(n_cols - 1)].reshape(
            (n_rows - 1, 2 * n_cols - 1)
        )
        tail_nodes[:, : n_cols - 1] = nodes[:-1, 1:]
        tail_nodes[:, n_cols - 1 :] = nodes[1:, :]

        return nodes_at_link.T

    @staticmethod
    def horizontal_links(shape):
        n_rows, n_cols = shape
        horizontal_links = np.empty((n_rows, n_cols - 1), dtype=int)
        horizontal_links[:, :] = np.arange(n_cols - 1)
        horizontal_links[:, :] += np.arange(n_rows).reshape((n_rows, 1)) * (
            2 * n_cols - 1
        )

        return horizontal_links.reshape(-1)

    @staticmethod
    def vertical_links(shape):
        n_rows, n_cols = shape

        vertical_links = np.empty((n_rows - 1, n_cols), dtype=int)
        vertical_links[:, :] = np.arange(n_cols) + n_cols - 1
        vertical_links[:, :] += np.arange(n_rows - 1).reshape((n_rows - 1, 1)) * (
            2 * n_cols - 1
        )

        return vertical_links.reshape(-1)

    @staticmethod
    def perimeter_nodes(shape):
        n_rows, n_cols = shape
        (
            northeast,
            northwest,
            southwest,
            southeast,
        ) = StructuredQuadLayout.corner_nodes(shape)

        return np.concatenate(
            (
                np.arange(southeast, northeast, n_cols),
                np.arange(northeast, northwest, -1),
                np.arange(northwest, southwest, -n_cols),
                np.arange(southwest, southeast, 1),
            )
        )

    @staticmethod
    def links_at_node(shape):
        n_rows, n_cols = shape

        links_at_node = np.empty((n_rows * n_cols, 4), dtype=int)

        east_links_at_node = links_at_node[:, 0].reshape((n_rows, n_cols))[:, :-1]
        east_links_at_node[:] = StructuredQuadLayoutPython.horizontal_links(
            shape
        ).reshape((n_rows, n_cols - 1))
        west_links_at_node = links_at_node[:, 2].reshape((n_rows, n_cols))[:, 1:]
        west_links_at_node[:] = StructuredQuadLayoutPython.horizontal_links(
            shape
        ).reshape((n_rows, n_cols - 1))

        north_links_at_node = links_at_node[:, 1].reshape((n_rows, n_cols))[:-1, :]
        north_links_at_node[:] = StructuredQuadLayoutPython.vertical_links(
            shape
        ).reshape((n_rows - 1, n_cols))
        south_links_at_node = links_at_node[:, 3].reshape((n_rows, n_cols))[1:, :]
        south_links_at_node[:] = StructuredQuadLayoutPython.vertical_links(
            shape
        ).reshape((n_rows - 1, n_cols))

        (
            northeast,
            northwest,
            southwest,
            southeast,
        ) = StructuredQuadLayout.corner_nodes(shape)
        bottom = slice(southwest, southeast + 1)
        top = slice(northwest, northeast + 1)
        left = slice(southwest, northwest + n_cols, n_cols)
        right = slice(southeast, northeast + n_cols, n_cols)

        for col, edge in enumerate((right, top, left, bottom)):
            links_at_node[edge, col] = -1

        return links_at_node

    @staticmethod
    def patches_at_link(shape):
        n_rows, n_cols = shape
        n_links = shape[0] * (shape[1] - 1) + (shape[0] - 1) * shape[1]
        n_patches = (n_rows - 1) * (n_cols - 1)
        patches = np.arange(n_patches, dtype=int).reshape((n_rows - 1, n_cols - 1))

        patches_at_link = np.empty((2, n_links), dtype=int)
        patches_at_link[0, : n_cols - 1] = -1
        patches_at_link[1, -(n_cols - 1) :] = -1

        right = (0, slice(n_cols - 1, n_links))
        left = (1, slice(0, -(n_cols - 1)))
        for edge in (right, left):
            edge_patches = patches_at_link[edge].reshape((n_rows - 1, 2 * n_cols - 1))
            edge_patches[:, : n_cols - 1] = patches
            edge_patches[:, n_cols - 1] = -1
            edge_patches[:, -(n_cols - 1) :] = patches

        return patches_at_link.T

    @staticmethod
    def link_dirs_at_node(shape):
        n_rows, n_cols = shape
        (
            northeast,
            northwest,
            southwest,
            southeast,
        ) = StructuredQuadLayout.corner_nodes(shape)

        link_dirs_at_node = np.empty((n_rows * n_cols, 4), dtype=np.int8)
        link_dirs_at_node[:, 0] = -1
        link_dirs_at_node[:, 1] = -1
        link_dirs_at_node[:, 2] = 1
        link_dirs_at_node[:, 3] = 1

        bottom = slice(southwest, southeast + 1)
        top = slice(northwest, northeast + 1)
        left = slice(southwest, northwest + n_cols, n_cols)
        right = slice(southeast, northeast + n_cols, n_cols)

        for col, edge in enumerate((right, top, left, bottom)):
            link_dirs_at_node[edge, col] = 0

        return link_dirs_at_node

    @staticmethod
    def patches_at_node(shape):
        n_rows, n_cols = shape

        patches_at_node = np.empty((4, n_rows * n_cols), dtype=int)

        ne = (slice(n_rows - 1), slice(n_cols - 1))
        nw = (slice(n_rows - 1), slice(1, n_cols))
        sw = (slice(1, n_rows), slice(1, n_cols))
        se = (slice(1, n_rows), slice(n_cols - 1))

        patches = np.arange((n_rows - 1) * (n_cols - 1), dtype=int).reshape(
            (n_rows - 1, n_cols - 1)
        )
        for col, nodes in enumerate((ne, nw, sw, se)):
            patches_at_node[col].reshape((n_rows, n_cols))[nodes] = patches

        (
            northeast,
            northwest,
            southwest,
            southeast,
        ) = StructuredQuadLayout.corner_nodes(shape)
        bottom = slice(southwest, southeast + 1)
        top = slice(northwest, northeast + 1)
        left = slice(southwest, northwest + n_cols, n_cols)
        right = slice(southeast, northeast + n_cols, n_cols)

        patches_at_node[0, right] = -1
        patches_at_node[0, top] = -1
        patches_at_node[1, left] = -1
        patches_at_node[1, top] = -1
        patches_at_node[2, left] = -1
        patches_at_node[2, bottom] = -1
        patches_at_node[3, right] = -1
        patches_at_node[3, bottom] = -1

        return patches_at_node.T


class StructuredQuadGraphTopology:
    _layout = StructuredQuadLayoutCython

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
    @read_only_array
    def nodes(self):
        """A shaped array of node ids.

        Returns
        -------
        ndarray
            Node IDs in an array shaped as *number_of_node_rows* by
            *number_of_node_columns*.
        """
        return np.arange(self.shape[0] * self.shape[1]).reshape(self.shape)

    @property
    @lru_cache()
    @read_only_array
    def nodes_at_right_edge(self):
        return np.arange(self.shape[1] - 1, np.prod(self.shape), self.shape[1])

    @property
    @lru_cache()
    @read_only_array
    def nodes_at_top_edge(self):
        return np.arange(self.number_of_nodes - self.shape[1], np.prod(self.shape))

    @property
    @lru_cache()
    @read_only_array
    def nodes_at_left_edge(self):
        return np.arange(0, np.prod(self.shape), self.shape[1])

    @property
    @lru_cache()
    @read_only_array
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

    @property
    @lru_cache()
    @read_only_array
    def nodes_at_link(self):
        return self._layout.nodes_at_link(self.shape)

    @property
    @lru_cache()
    def horizontal_links(self):
        return self._layout.horizontal_links(self.shape)

    @property
    @lru_cache()
    def vertical_links(self):
        return self._layout.vertical_links(self.shape)

    @property
    def corner_nodes(self):
        n_rows, n_cols = self.shape
        return np.asarray(
            (n_rows * n_cols - 1, (n_rows - 1) * n_cols, 0, n_cols - 1), dtype=int
        )

    @property
    @lru_cache()
    def perimeter_nodes(self):
        return self._layout.perimeter_nodes(self.shape)

    @property
    @lru_cache()
    def links_at_node(self):
        return self._layout.links_at_node(self.shape)

    @property
    @lru_cache()
    def link_dirs_at_node(self):
        return self._layout.link_dirs_at_node(self.shape)

    @property
    @lru_cache()
    @read_only_array
    def patches_at_link(self):
        return self._layout.patches_at_link(self.shape)

    @property
    @lru_cache()
    @read_only_array
    def patches_at_node(self):
        return self._layout.patches_at_node(self.shape)


class StructuredQuadGraphExtras(StructuredQuadGraphTopology, Graph):
    def __init__(self, node_y_and_x, sort=False):
        StructuredQuadGraphTopology.__init__(self, node_y_and_x[0].shape)
        Graph.__init__(
            self,
            node_y_and_x,
            links=StructuredQuadLayoutCython.nodes_at_link(self.shape),
            patches=StructuredQuadLayoutCython.links_at_patch(self.shape),
            sort=sort,
        )

    @property
    def nodes_at_link(self):
        return self.ds["nodes_at_link"].values


class StructuredQuadGraph(StructuredQuadGraphExtras):
    def __init__(self, coords, shape=None, sort=False):
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

        StructuredQuadGraphExtras.__init__(self, (node_y, node_x), sort=sort)

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

    def __init__(self, nodes, sort=False):
        rows = np.asarray(nodes[0], dtype=float)
        cols = np.asarray(nodes[1], dtype=float)
        node_y_and_x = np.meshgrid(rows, cols, indexing="ij")

        StructuredQuadGraphExtras.__init__(self, node_y_and_x, sort=sort)

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

    def __init__(self, shape, spacing=1.0, origin=0.0, sort=False):
        spacing = np.broadcast_to(spacing, 2)
        origin = np.broadcast_to(origin, 2)

        rows = np.arange(shape[0], dtype=float) * spacing[0] + origin[0]
        cols = np.arange(shape[1], dtype=float) * spacing[1] + origin[1]

        node_y_and_x = np.meshgrid(rows, cols, indexing="ij")

        StructuredQuadGraphExtras.__init__(self, node_y_and_x, sort=sort)

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

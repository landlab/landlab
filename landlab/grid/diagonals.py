#! /usr/bin/env python
import numpy as np

from ..utils.decorators import cache_result_in_object, make_return_array_immutable
from .decorators import return_readonly_id_array
from .linkstatus import ACTIVE_LINK, set_status_at_link


def create_nodes_at_diagonal(shape, out=None):
    """Create array of tail and head nodes for diagonals.

    Parameters
    ----------
    shape : tuple of *(n_rows, n_cols)*
        Shape as number of node rows and node columns.
    out : ndarray of shape *(n_diagonals, 2)*, optional
        Output buffer to place nodes at each diagonal.

    Returns
    -------
    out : ndarray of shape *(n_diagonals, 2)*
        Tail and head node for each diagonal.

    Examples
    --------
    >>> from landlab.grid.diagonals import create_nodes_at_diagonal
    >>> create_nodes_at_diagonal((3, 4)) # doctest: +NORMALIZE_WHITESPACE
    array([[ 0, 5], [ 1, 4], [ 1,  6], [ 2, 5], [ 2,  7], [ 3,  6],
           [ 4, 9], [ 5, 8], [ 5, 10], [ 6, 9], [ 6, 11], [ 7, 10]])
    """
    shape = np.asarray(shape)
    n_diagonals = np.prod(shape - 1) * 2
    n_nodes = np.prod(shape)
    if out is None:
        out = np.empty((n_diagonals, 2), dtype=int)

    nodes = np.arange(n_nodes).reshape(shape)

    out[::2, 0] = nodes[:-1, :-1].flat
    out[::2, 1] = nodes[1:, 1:].flat
    out[1::2, 0] = nodes[:-1, 1:].flat
    out[1::2, 1] = nodes[1:, :-1].flat

    return out


def create_diagonals_at_node(shape, out=None):
    """Create array of diagonals at node.

    Create an array that gives the diagonal touching each node
    for a structured grid of quadrilaterals. For each node, links
    are ordered clockwise from axis=1.

    Parameters
    ----------
    shape : tuple of *(n_rows, n_cols)*
        Shape as number of node rows and node columns.
    out : ndarray of shape *(n_nodes, 4)*, optional
        Output buffer to place diagonal ids at each node.

    Returns
    -------
    out : ndarray of shape *(n_nodes, 4)*
        Diagonals at node with -1 for missing diagonals.

    Examples
    --------
    >>> from landlab.grid.diagonals import create_diagonals_at_node
    >>> create_diagonals_at_node((3, 4))
    array([[ 0, -1, -1, -1],
           [ 2,  1, -1, -1],
           [ 4,  3, -1, -1],
           [-1,  5, -1, -1],
           [ 6, -1, -1,  1],
           [ 8,  7,  0,  3],
           [10,  9,  2,  5],
           [-1, 11,  4, -1],
           [-1, -1, -1,  7],
           [-1, -1,  6,  9],
           [-1, -1,  8, 11],
           [-1, -1, 10, -1]])
    """
    shape = np.asarray(shape)
    n_diagonals = np.prod(shape - 1) * 2
    n_nodes = np.prod(shape)
    if out is None:
        out = np.full((n_nodes, 4), -1, dtype=int)

    diagonals = np.full(shape + 1, -1, dtype=int)

    diagonals[1:-1, 1:-1] = np.arange(0, n_diagonals, 2).reshape(shape - 1)
    out[:, 0] = diagonals[1:, 1:].flat
    out[:, 2] = diagonals[:-1, :-1].flat

    diagonals[1:-1, 1:-1] = np.arange(1, n_diagonals, 2).reshape(shape - 1)
    out[:, 1] = diagonals[1:, :-1].flat
    out[:, 3] = diagonals[:-1, 1:].flat

    return out


class DiagonalsMixIn(object):

    """Add diagonals to a structured quad grid."""

    @property
    @cache_result_in_object()
    def number_of_diagonals(self):
        """Number of diagonals in the grid.

        Examples
        --------
        >>> from landlab import RasterModelGrid
        >>> grid = RasterModelGrid((4, 3))
        >>> grid.number_of_diagonals
        12
        """
        return 2 * np.prod(np.asarray(self.shape) - 1)

    @property
    @cache_result_in_object()
    @make_return_array_immutable
    def diagonals_at_node(self):
        """Diagonals attached to nodes.

        Returns
        -------
        ndarray of int, shape `(n_nodes, 4)`
            Diagonals at each node.

        Examples
        --------
        >>> from landlab import RasterModelGrid
        >>> grid = RasterModelGrid((4, 3))
        >>> grid.diagonals_at_node.shape == (grid.number_of_nodes, 4)
        True
        >>> grid.diagonals_at_node
        array([[ 0, -1, -1, -1], [ 2,  1, -1, -1], [-1,  3, -1, -1],
               [ 4, -1, -1,  1], [ 6,  5,  0,  3], [-1,  7,  2, -1],
               [ 8, -1, -1,  5], [10,  9,  4,  7], [-1, 11,  6, -1],
               [-1, -1, -1,  9], [-1, -1,  8, 11], [-1, -1, 10, -1]])

        LLCATS: NINF LINF CONN
        """
        return create_diagonals_at_node(self.shape)

    @property
    @cache_result_in_object()
    @make_return_array_immutable
    def diagonal_dirs_at_node(self):
        """Directions of diagonals attached to nodes.

        Array of diagonal directions relative to each node. If the
        diagonal is directed away from the node the direction is -1, if the
        diagonal is directed toward the node its direction is 1. Each
        node is assumed to have exactly four diagonals attached to it.
        However, perimeter nodes will have few diagonals (corner ndoes
        will only have one diagonal and edge nodes two). In such cases,
        placeholders of 0 are used.

        Examples
        --------
        >>> from landlab import RasterModelGrid
        >>> grid = RasterModelGrid((4, 3))
        >>> grid.diagonal_dirs_at_node.shape == (grid.number_of_nodes, 4)
        True
        >>> grid.diagonal_dirs_at_node
        array([[-1,  0,  0,  0], [-1, -1,  0,  0], [ 0, -1,  0,  0],
               [-1,  0,  0,  1], [-1, -1,  1,  1], [ 0, -1,  1,  0],
               [-1,  0,  0,  1], [-1, -1,  1,  1], [ 0, -1,  1,  0],
               [ 0,  0,  0,  1], [ 0,  0,  1,  1], [ 0,  0,  1,  0]],
               dtype=int8)

        The lower-right corner node only has one diagonal.

        >>> grid.diagonal_dirs_at_node[2]
        array([ 0, -1,  0,  0], dtype=int8)

        A node on one of the edges has two diagonals.

        >>> grid.diagonal_dirs_at_node[3]
        array([-1,  0,  0,  1], dtype=int8)
        """
        diagonals_at_node = self.diagonals_at_node
        dirs_at_node = np.zeros_like(diagonals_at_node, dtype=np.int8)
        dirs_at_node[diagonals_at_node >= 0] = 1
        dirs_at_node[:, (0, 1)] *= -1

        return dirs_at_node

    @property
    @cache_result_in_object()
    @make_return_array_immutable
    def diagonal_adjacent_nodes_at_node(self):
        """Get adjacent nodes along diagonals.

        Order is the landlab standard, counterclockwise starting from
        east.

        Examples
        --------
        >>> from landlab import RasterModelGrid
        >>> grid = RasterModelGrid((4, 3))
        >>> grid.diagonal_adjacent_nodes_at_node
        array([[ 4, -1, -1, -1], [ 5,  3, -1, -1], [-1,  4, -1, -1],
               [ 7, -1, -1,  1], [ 8,  6,  0,  2], [-1,  7,  1, -1],
               [10, -1, -1,  4], [11,  9,  3,  5], [-1, 10,  4, -1],
               [-1, -1, -1,  7], [-1, -1,  6,  8], [-1, -1,  7, -1]])

        LLCATS: NINF CONN
        """
        node_is_at_tail = np.choose(
            self.diagonal_dirs_at_node + 1, np.array((1, -1, 0), dtype=np.int8)
        )
        out = self.nodes_at_diagonal[self.diagonals_at_node, node_is_at_tail]
        out[node_is_at_tail == -1] = -1

        return out

    @property
    @cache_result_in_object()
    @make_return_array_immutable
    def d8_adjacent_nodes_at_node(self):
        return np.vstack(
            (
                super(DiagonalsMixIn, self).adjacent_nodes_at_node,
                self.diagonal_adjacent_nodes_at_node,
            )
        )

    @property
    @cache_result_in_object()
    @make_return_array_immutable
    def nodes_at_diagonal(self):
        """Nodes at diagonal tail and head.

        Examples
        --------
        >>> from landlab import RasterModelGrid
        >>> grid = RasterModelGrid((4, 3))
        >>> grid.nodes_at_diagonal.shape == (grid.number_of_diagonals, 2)
        True
        >>> grid.nodes_at_diagonal
        array([[ 0,  4], [ 1,  3], [ 1,  5], [ 2,  4],
               [ 3,  7], [ 4,  6], [ 4,  8], [ 5,  7],
               [ 6, 10], [ 7,  9], [ 7, 11], [ 8, 10]])

        The first column is diagonal tail and the second the head.

        >>> grid.diagonals_at_node[3]
        array([ 4, -1, -1,  1])
        >>> grid.nodes_at_diagonal[(4, 1), ]
        array([[3, 7],
               [1, 3]])
        >>> grid.diagonal_dirs_at_node[3]
        array([-1,  0,  0,  1], dtype=int8)
        """
        return create_nodes_at_diagonal(self.shape)

    @property
    @cache_result_in_object()
    def number_of_d8(self):
        """Number of links and diagonals.

        Examples
        --------
        >>> from landlab import RasterModelGrid
        >>> grid = RasterModelGrid((4, 3))
        >>> grid.number_of_d8
        29
        """
        return super(DiagonalsMixIn, self).number_of_links + self.number_of_diagonals

    @property
    @cache_result_in_object()
    @make_return_array_immutable
    def nodes_at_d8(self):
        return np.vstack((self.nodes_at_link, self.nodes_at_diagonal))

    @property
    @cache_result_in_object()
    @make_return_array_immutable
    def d8s_at_node(self):
        """Links and diagonals attached to nodes.

        Returns
        -------
        ndarray of int, shape `(n_nodes, 8)`
            Links and diagonals at each node.

        Examples
        --------
        >>> import numpy as np
        >>> from landlab import RasterModelGrid

        >>> grid = RasterModelGrid((3, 4))
        >>> grid.d8s_at_node.shape == (grid.number_of_nodes, 8)
        True
        >>> grid.d8s_at_node
        array([[ 0,  3, -1, -1, 17, -1, -1, -1],
               [ 1,  4,  0, -1, 19, 18, -1, -1],
               [ 2,  5,  1, -1, 21, 20, -1, -1],
               [-1,  6,  2, -1, -1, 22, -1, -1],
               [ 7, 10, -1,  3, 23, -1, -1, 18],
               [ 8, 11,  7,  4, 25, 24, 17, 20],
               [ 9, 12,  8,  5, 27, 26, 19, 22],
               [-1, 13,  9,  6, -1, 28, 21, -1],
               [14, -1, -1, 10, -1, -1, -1, 24],
               [15, -1, 14, 11, -1, -1, 23, 26],
               [16, -1, 15, 12, -1, -1, 25, 28],
               [-1, -1, 16, 13, -1, -1, 27, -1]])
        >>> np.all(grid.d8s_at_node[:, :4] == grid.links_at_node)
        True

        >>> diagonals_at_node = grid.d8s_at_node[:, 4:] - grid.number_of_links
        >>> diagonals_at_node[grid.d8s_at_node[:, 4:] == -1] = -1
        >>> np.all(diagonals_at_node == grid.diagonals_at_node)
        True

        LLCATS: NINF LINF CONN
        """
        diagonals_at_node = self.diagonals_at_node.copy()
        diagonals_at_node[diagonals_at_node >= 0] += self.number_of_links
        return np.hstack((super(DiagonalsMixIn, self).links_at_node, diagonals_at_node))

    @property
    @cache_result_in_object()
    @make_return_array_immutable
    def d8_dirs_at_node(self):
        return np.hstack(
            (super(DiagonalsMixIn, self).link_dirs_at_node, self.diagonal_dirs_at_node)
        )

    @property
    # @cache_result_in_object()
    @make_return_array_immutable
    def d8_status_at_node(self):
        return self.status_at_d8[self.d8s_at_node]

    @property
    @cache_result_in_object()
    @make_return_array_immutable
    def length_of_diagonal(self):
        return np.sqrt(
            np.power(np.diff(self.xy_of_node[self.nodes_at_diagonal], axis=1), 2.0).sum(
                axis=2
            )
        ).flatten()

    @property
    @cache_result_in_object()
    @make_return_array_immutable
    def length_of_d8(self):
        """Length of links and diagonals.

        Return the lengths links and diagonals in the grid. Links are
        listed first and then diagonals.

        Returns
        -------
        ndarray of float
            Link lengths.

        Examples
        --------
        >>> from landlab import RasterModelGrid
        >>> grid = RasterModelGrid((3, 3), xy_spacing=(4, 3))

        >>> grid.length_of_link
        array([ 4.,  4.,  3.,  3.,  3.,  4.,  4.,  3.,  3.,  3.,  4.,  4.])

        >>> grid.length_of_d8 # doctest: +NORMALIZE_WHITESPACE
        array([ 4.,  4.,  3.,  3.,  3.,
                4.,  4.,  3.,  3.,  3.,
                4.,  4.,  5.,  5.,  5.,
                5.,  5.,  5.,  5.,  5.])

        LLCATS: LINF MEAS
        """
        return np.hstack(
            (super(DiagonalsMixIn, self).length_of_link, self.length_of_diagonal)
        )

    def reset_status_at_node(self):
        super(DiagonalsMixIn, self).reset_status_at_node()
        attrs = [
            "_status_at_diagonal",
            "_diagonal_status_at_node",
            "_active_diagonals",
            "_active_diagonal_dirs_at_node",
            "_status_at_d8",
            "_active_d8",
            "_active_d8_dirs_at_node",
        ]

        for attr in attrs:
            try:
                del self.__dict__[attr]
            except KeyError:
                pass

    @property
    @cache_result_in_object()
    @make_return_array_immutable
    def status_at_diagonal(self):
        """Status at diagonals.

        Examples
        --------
        >>> from landlab import RasterModelGrid
        >>> grid = RasterModelGrid((4, 3))

        An inactive link (or diagonal) is one that joins two
        boundary nodes or has one end that is a closed boundary.

        >>> grid.status_at_node
        array([1, 1, 1,
               1, 0, 1,
               1, 0, 1,
               1, 1, 1], dtype=uint8)
        >>> grid.BC_NODE_IS_CORE, grid.BC_NODE_IS_FIXED_VALUE
        (0, 1)

        Diagonals that connect two boundary nodes are inactive.

        >>> grid.status_at_diagonal
        array([0, 4, 4, 0,
               0, 0, 0, 0,
               4, 0, 0, 4], dtype=uint8)
        >>> grid.BC_LINK_IS_ACTIVE, grid.BC_LINK_IS_INACTIVE
        (0, 4)

        By setting a node to closed that wasn't before, a new link
        become inactive.

        >>> grid.status_at_node[5] = grid.BC_NODE_IS_CLOSED
        >>> grid.status_at_diagonal
        array([0, 4, 4, 0,
               0, 0, 0, 4,
               4, 0, 0, 4], dtype=uint8)
        """
        return set_status_at_link(self.status_at_node[self.nodes_at_diagonal])

    @property
    @cache_result_in_object()
    @make_return_array_immutable
    def diagonal_status_at_node(self):
        return self.status_at_diagonal[self.diagonals_at_node]

    @property
    @cache_result_in_object()
    @return_readonly_id_array
    def active_diagonals(self):
        return np.where(self.status_at_diagonal == ACTIVE_LINK)[0]

    @property
    @cache_result_in_object()
    @make_return_array_immutable
    def active_diagonal_dirs_at_node(self):
        return np.choose(
            self.diagonal_status_at_node == ACTIVE_LINK, (0, self.diagonal_dirs_at_node)
        )

    @property
    @cache_result_in_object()
    @make_return_array_immutable
    def status_at_d8(self):
        return np.hstack(
            (super(DiagonalsMixIn, self).status_at_link, self.status_at_diagonal)
        )

    @property
    @cache_result_in_object()
    @return_readonly_id_array
    def active_d8(self):
        return np.where(self.status_at_d8 == ACTIVE_LINK)[0]

    @property
    @cache_result_in_object()
    @make_return_array_immutable
    def active_d8_dirs_at_node(self):
        return np.choose(
            self.d8_status_at_node == ACTIVE_LINK, (0, self.d8_dirs_at_node)
        )

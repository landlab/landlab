#! /usr/bin/env python
import numpy as np


from .base import (ACTIVE_LINK, INACTIVE_LINK, FIXED_LINK,
                   CLOSED_BOUNDARY, CORE_NODE, FIXED_GRADIENT_BOUNDARY,
                   FIXED_VALUE_BOUNDARY)
from ..utils.decorators import (cache_result_in_object,
                                make_return_array_immutable)


def is_fixed_link(node_status_at_link):
    """Find links that are fixed.

    A link is fixed if it connects a core node with a fixed value
    boundary node.

    Parameters
    ----------
    node_status_at_link : ndarray of int, shape `(n_links, 2)`
        Node status a link tail and head.

    Returns
    -------
    ndarray of bool, shape `(n_links, )`
        True if link is fixed.

    Examples
    --------
    >>> from landlab.grid.diagonals import is_fixed_link
    >>> from landlab import CORE_NODE, FIXED_GRADIENT_BOUNDARY
    >>> is_fixed_link([CORE_NODE, FIXED_GRADIENT_BOUNDARY])
    array([ True], dtype=bool)

    >>> from landlab import FIXED_VALUE_BOUNDARY
    >>> is_fixed_link([CORE_NODE, FIXED_VALUE_BOUNDARY])
    array([False], dtype=bool)

    >>> is_fixed_link([[FIXED_GRADIENT_BOUNDARY, CORE_NODE],
    ...                [CORE_NODE, CORE_NODE]])
    array([ True, False], dtype=bool)
    """
    node_status_at_link = np.asarray(node_status_at_link).reshape((-1, 2))

    is_core_node = node_status_at_link == CORE_NODE
    is_fixed_gradient_node = node_status_at_link == FIXED_GRADIENT_BOUNDARY

    return ((is_core_node[:, 0] & is_fixed_gradient_node[:, 1]) |
            (is_fixed_gradient_node[:, 0] & is_core_node[:, 1]))


def is_inactive_link(node_status_at_link):
    """Find links that are inactive.

    A link is inactive if it connects two boundary nodes or one of
    its nodes is closed.

    Parameters
    ----------
    node_status_at_link : ndarray of int, shape `(n_links, 2)`
        Node status a link tail and head.

    Returns
    -------
    ndarray of bool, shape `(n_links, )`
        True if link is isactive.

    Examples
    --------
    >>> from landlab.grid.diagonals import is_inactive_link
    >>> from landlab import CORE_NODE, FIXED_GRADIENT_BOUNDARY
    >>> is_inactive_link([CORE_NODE, CLOSED_BOUNDARY])
    array([ True], dtype=bool)

    >>> from landlab import FIXED_VALUE_BOUNDARY
    >>> is_inactive_link([FIXED_GRADIENT_BOUNDARY, FIXED_VALUE_BOUNDARY])
    array([ True], dtype=bool)

    >>> is_inactive_link([[FIXED_GRADIENT_BOUNDARY, CLOSED_BOUNDARY],
    ...                   [CORE_NODE, CORE_NODE]])
    array([ True, False], dtype=bool)
    """
    node_status_at_link = np.asarray(node_status_at_link).reshape((-1, 2))

    is_core = node_status_at_link == CORE_NODE
    is_fixed_value = node_status_at_link == FIXED_VALUE_BOUNDARY
    is_fixed_gradient = node_status_at_link == FIXED_GRADIENT_BOUNDARY
    is_closed = node_status_at_link == CLOSED_BOUNDARY
    is_boundary_node = is_fixed_value | is_fixed_gradient | is_closed

    return ((is_boundary_node[:, 0] & is_boundary_node[:, 1]) |
            (is_closed[:, 0] & is_core[:, 1]) |
            (is_core[:, 0] & is_closed[:, 1]))


def is_active_link(node_status_at_link):
    """Find links that are active.

    A link is active if it connects a core node with another core
    node or a fixed value boundary.

    Parameters
    ----------
    node_status_at_link : ndarray of int, shape `(n_links, 2)`
        Node status a link tail and head.

    Returns
    -------
    ndarray of bool, shape `(n_links, )`
        True if link is isactive.

    Examples
    --------
    >>> from landlab.grid.diagonals import is_active_link
    >>> from landlab import CORE_NODE, FIXED_GRADIENT_BOUNDARY
    >>> is_active_link([CORE_NODE, FIXED_GRADIENT_BOUNDARY])
    array([False], dtype=bool)

    >>> from landlab import FIXED_VALUE_BOUNDARY
    >>> is_active_link([CORE_NODE, FIXED_VALUE_BOUNDARY])
    array([ True], dtype=bool)

    >>> is_active_link([[FIXED_GRADIENT_BOUNDARY, CORE_NODE],
    ...                 [CORE_NODE, CORE_NODE]])
    array([False, True], dtype=bool)
    """
    node_status_at_link = np.asarray(node_status_at_link).reshape((-1, 2))

    is_core_node = node_status_at_link == CORE_NODE
    is_fixed_value_node = node_status_at_link == FIXED_VALUE_BOUNDARY
    return (
        (is_core_node[:, 0] & is_core_node[:, 1]) |
        (is_core_node[:, 0] & is_fixed_value_node[:, 1]) |
        (is_fixed_value_node[:, 0] & is_core_node[:, 1])
    )


def set_status_at_link(node_status_at_link, out=None):
    n_links = len(node_status_at_link)

    if out is None:
        out = np.full(n_links, 255, dtype=np.uint8)

    _is_fixed_link = is_fixed_link(node_status_at_link)
    _is_active_link = is_active_link(node_status_at_link)
    _is_inactive_link = is_inactive_link(node_status_at_link)

    assert np.all(np.sum(np.vstack((_is_active_link, _is_inactive_link,
                                    _is_fixed_link)), axis=0) == 1)

    out[_is_inactive_link] = INACTIVE_LINK
    out[_is_active_link] = ACTIVE_LINK
    out[_is_fixed_link] = FIXED_LINK

    return out


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


def create_diagonal_dirs_at_node(shape, out=None):
    """Create array of diagonals directions at node.

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
        out = np.full((n_nodes, 4), -1, dtype=int8)

    dirs = np.zeros(shape + 1, dtype=int8)

    dirs[1:-1, 1:-1] = np.arange(0, n_diagonals, 2).reshape(shape - 1)
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
    def diagonal_dirs_at_node(self):
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

        LLCATS: DEPR NINF CONN
        """
        node_is_at_tail = np.choose(self.diagonal_dirs_at_node + 1,
                                    np.array((1, -1, 0), dtype=np.int8))
        out = self.nodes_at_diagonal[self.diagonals_at_node, node_is_at_tail]
        out[node_is_at_tail == -1] = -1

        return out

    @property
    @cache_result_in_object()
    def d8_adjacent_nodes_at_node(self):
        return np.vstack((super(DiagonalsMixIn, self).adjacent_nodes_at_node,
                          self.diagonal_adjacent_nodes_at_node))

    @property
    @cache_result_in_object()
    def diagonal_status_at_node(self):
        return self.status_at_diagonal[self.diagonals_at_node]

    @property
    @cache_result_in_object()
    def nodes_at_diagonal(self):
        return create_nodes_at_diagonal(self.shape)

    @property
    @cache_result_in_object(cache_as='_status_at_diagonal')
    def status_at_diagonal(self):
        return set_status_at_link(self.status_at_node[self.nodes_at_diagonal])

    @property
    @cache_result_in_object()
    def number_of_d8(self):
        return (super(DiagonalsMixIn, self).number_of_links +
                self.number_of_diagonals)

    @property
    @cache_result_in_object()
    def nodes_at_d8(self):
        return np.vstack((self.nodes_at_link, self.nodes_at_diagonal))

    @property
    @cache_result_in_object()
    def status_at_d8(self):
        return np.hstack((super(DiagonalsMixIn, self).status_at_link,
                          self.status_at_diagonal))

    @property
    @cache_result_in_object()
    @make_return_array_immutable
    def d8_at_node(self):
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
        >>> grid.d8_at_node.shape == (grid.number_of_nodes, 8)
        True
        >>> grid.d8_at_node
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
        >>> np.all(grid.d8_at_node[:, :4] == grid.links_at_node)
        True

        >>> diagonals_at_node = grid.d8_at_node[:, 4:] - grid.number_of_links
        >>> diagonals_at_node[grid.d8_at_node[:, 4:] == -1] = -1
        >>> np.all(diagonals_at_node == grid.diagonals_at_node)
        True

        LLCATS: NINF LINF CONN
        """
        diagonals_at_node = self.diagonals_at_node.copy()
        diagonals_at_node[diagonals_at_node >= 0] += self.number_of_links
        return np.hstack((super(DiagonalsMixIn, self).links_at_node,
                          diagonals_at_node))
                          # self.diagonals_at_node + self.number_of_links))

    @property
    @cache_result_in_object()
    def d8_dirs_at_node(self):
        return np.hstack((super(DiagonalsMixIn, self).link_dirs_at_node,
                          self.diagonal_dirs_at_node))

    @property
    @cache_result_in_object()
    def length_of_diagonal(self):
        return np.sqrt(
            np.power(np.diff(self.xy_of_node[self.nodes_at_diagonal], axis=1),
                     2.).sum(axis=2)).flatten()

    @property
    @cache_result_in_object()
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
        >>> grid = RasterModelGrid((3, 3), spacing=(3, 4))

        >>> grid.length_of_link
        array([ 4.,  4.,  3.,  3.,  3.,  4.,  4.,  3.,  3.,  3.,  4.,  4.])

        >>> grid.length_of_d8 # doctest: +NORMALIZE_WHITESPACE
        array([ 4.,  4.,  3.,  3.,  3.,
                4.,  4.,  3.,  3.,  3.,
                4.,  4.,  5.,  5.,  5.,
                5.,  5.,  5.,  5.,  5.])

        LLCATS: LINF MEAS
        """
        return np.hstack((super(DiagonalsMixIn, self).length_of_link,
                          self.length_of_diagonal))


# from .raster import RasterModelGrid


# class DiagonalRasterModelGrid(DiagonalsMixIn, RasterModelGrid):

#     pass

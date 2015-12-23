import numpy as np
from . import nodes
from ..base import CORE_NODE, FIXED_GRADIENT_BOUNDARY, FIXED_VALUE_BOUNDARY
from ..unstructured.links import LinkGrid
from ...core.utils import as_id_array


def neighbors_at_link(shape, links):
    """Get neighbor links.

    Examples
    --------
    >>> import numpy as np
    >>> from landlab.grid.structured_quad.links import neighbors_at_link

    >>> neighbors_at_link((3, 2), np.arange(7)) # doctest: +NORMALIZE_WHITESPACE
    array([[-1,  3, -1, -1],
           [ 2,  4, -1, -1], [-1,  5,  1, -1],
           [-1,  6, -1,  0],
           [ 5,  7, -1,  1], [-1, -1,  4,  2],
           [-1, -1, -1,  3]])
    """
    from .cfuncs import _neighbors_at_link

    links = np.asarray(links, dtype=int)
    out = np.full((links.size, 4), -1, dtype=int)
    _neighbors_at_link(links, shape, out)
    return out


def shape_of_vertical_links(shape):
    """Shape of vertical link grid.

    Number of rows and columns of *vertical* links that connect nodes in a
    structured grid of quadrilaterals.

    Parameters
    ----------
    shape : tuple of int
        Shape of grid of nodes.

    Returns
    -------
    tuple of int :
        Shape of the vertical links in grid.

    Examples
    --------
    >>> from landlab.grid.structured_quad.links import shape_of_vertical_links
    >>> shape_of_vertical_links((3, 4))
    (2, 4)
    """
    return (shape[0] - 1, shape[1])


def shape_of_horizontal_links(shape):
    """Shape of horizontal link grid.

    Number of rows and columns of *horizontal* links that connect nodes in a
    structured grid of quadrilaterals.

    Parameters
    ----------
    shape : tuple of int
        Shape of grid of nodes.

    Returns
    -------
    tuple of int :
        Shape of the horizontal links in grid.

    Examples
    --------
    >>> from landlab.grid.structured_quad.links import (
    ...     shape_of_horizontal_links)
    >>> shape_of_horizontal_links((3, 4))
    (3, 3)
    """
    return (shape[0], shape[1] - 1)


def number_of_vertical_links(shape):
    """Number of vertical links in a structured quad grid.

    Parameters
    ----------
    shape : tuple of int
        Shape of grid of nodes.

    Returns
    -------
    int :
        Number of vertical links in grid.

    Examples
    --------
    >>> from landlab.grid.structured_quad.links import number_of_vertical_links
    >>> number_of_vertical_links((3, 4))
    8
    """
    return np.prod(shape_of_vertical_links(shape))


def number_of_horizontal_links(shape):
    """Number of horizontal links in a structured quad grid.

    Parameters
    ----------
    shape : tuple of int
        Shape of grid of nodes.

    Returns
    -------
    int :
        Number of horizontal links in grid.

    Examples
    --------
    >>> from landlab.grid.structured_quad.links import (
    ...     number_of_horizontal_links)
    >>> number_of_horizontal_links((3, 4))
    9
    """
    return np.prod(shape_of_horizontal_links(shape))


def number_of_links(shape):
    """Number of links in a structured quad grid.

    Parameters
    ----------
    shape : tuple of int
        Shape of grid of nodes.

    Returns
    -------
    int :
        Number of links in grid.

    Examples
    --------
    >>> from landlab.grid.structured_quad.links import number_of_links
    >>> number_of_links((3, 4))
    17
    """
    return number_of_vertical_links(shape) + number_of_horizontal_links(shape)


def vertical_link_ids(shape):
    """Vertical links in a structured quad grid.

    Parameters
    ----------
    shape : tuple of int
        Shape of grid of nodes.

    Returns
    -------
    (M, N) ndarray :
        Array of link IDs.

    Examples
    --------
    >>> from landlab.grid.structured_quad.links import vertical_link_ids
    >>> vertical_link_ids((3, 4))
    array([[ 3,  4,  5,  6],
           [10, 11, 12, 13]])
    """
    #link_ids = np.arange(number_of_vertical_links(shape), dtype=np.int)
    #return link_ids.reshape(shape_of_vertical_links(shape))
    a = shape[1] - 1  # num horiz links in each row
    num_links_per_row = 2*shape[1] - 1  # each row has C-1 horiz + C vert
    link_ids = np.zeros(shape_of_vertical_links(shape), dtype=np.int)
    for r in range(shape[0]-1):  # num rows - 1
        link_ids[r,:] = a + (r * num_links_per_row) + np.arange(shape[1])
    return link_ids


def horizontal_link_ids(shape):
    """Horizontal links in a structured quad grid.

    Parameters
    ----------
    shape : tuple of int
        Shape of grid of nodes.

    Returns
    -------
    (M, N) ndarray :
        Array of link IDs.

    Examples
    --------
    >>> from landlab.grid.structured_quad.links import horizontal_link_ids
    >>> horizontal_link_ids((3, 4))
    array([[ 0,  1,  2],
           [ 7,  8,  9],
           [14, 15, 16]])
    """
    num_links_per_row = 2*shape[1] - 1  # each row has C-1 horiz + C vert
    link_ids = np.zeros(shape_of_horizontal_links(shape), dtype=np.int)
    for r in range(shape[0]):  # number of rows
        link_ids[r,:] = (r * num_links_per_row) + np.arange(shape[1]-1)
    return link_ids


def number_of_links_per_node(shape):
    """Number of links touching each node.

    Parameters
    ----------
    shape : tuple of int
        Shape of grid of nodes.

    Returns
    -------
    ndarray :
        Array of number of links per node.

    Examples
    --------
    >>> from landlab.grid.structured_quad.links import (
    ...     number_of_links_per_node, number_of_in_links_per_node,
    ...     number_of_out_links_per_node)
    >>> number_of_links_per_node((3, 4))
    array([[2, 3, 3, 2],
           [3, 4, 4, 3],
           [2, 3, 3, 2]])
    >>> (number_of_in_links_per_node((3, 4)) +
    ...  number_of_out_links_per_node((3, 4)))
    array([[2, 3, 3, 2],
           [3, 4, 4, 3],
           [2, 3, 3, 2]])
    """
    link_count = np.empty(shape, np.int)
    link_count[1:-1, 1:-1] = 4
    link_count[(0, -1), 1:-1] = 3
    link_count[1:-1, (0, -1)] = 3
    link_count[(0, 0, -1, -1), (0, -1, 0, -1)] = 2
    return link_count


def number_of_in_links_per_node(shape):
    """Number of links entering each node.

    Parameters
    ----------
    shape : tuple of int
        Shape of grid of nodes.

    Returns
    -------
    ndarray :
        Array of number of in-links per node.

    Examples
    --------
    >>> from landlab.grid.structured_quad.links import (
    ...     number_of_in_links_per_node)
    >>> number_of_in_links_per_node((3, 4))
    array([[0, 1, 1, 1],
           [1, 2, 2, 2],
           [1, 2, 2, 2]])
    """
    link_count = np.empty(shape, np.int)
    link_count[1:, 1:] = 2
    link_count[0, 0] = 0
    link_count[0, 1:] = 1
    link_count[1:, 0] = 1
    return link_count


def number_of_out_links_per_node(shape):
    """Number of links leaving each node.

    Parameters
    ----------
    shape : tuple of int
        Shape of grid of nodes.

    Returns
    -------
    ndarray :
        Array of number of out-links per node.

    Examples
    --------
    >>> from landlab.grid.structured_quad.links import (
    ...     number_of_out_links_per_node)
    >>> number_of_out_links_per_node((3, 4))
    array([[2, 2, 2, 1],
           [2, 2, 2, 1],
           [1, 1, 1, 0]])
    """
    link_count = np.empty(shape, np.int)
    link_count[:-1, :-1] = 2
    link_count[-1, -1] = 0
    link_count[-1, :-1] = 1
    link_count[:-1, -1] = 1
    return link_count


def _node_out_link_ids(shape):
    """Links leaving each node.

    Parameters
    ----------
    shape : tuple of int
        Shape of grid of nodes.

    Returns
    -------
    tuple :
        Tuple of array of link IDs as (vertical_links, horizontal_links).

    Examples
    --------
    >>> from landlab.grid.structured_quad.links import _node_out_link_ids
    >>> (vert, horiz) = _node_out_link_ids((3, 4))
    >>> vert
    array([[ 3,  4,  5,  6],
           [10, 11, 12, 13],
           [-1, -1, -1, -1]])
    >>> horiz
    array([[ 0,  1,  2, -1],
           [ 7,  8,  9, -1],
           [14, 15, 16, -1]])
    """
    node_horizontal_link_ids = np.empty(shape, np.int)
    node_horizontal_link_ids[:, :-1] = horizontal_link_ids(shape)
    node_horizontal_link_ids[:, -1] = -1

    node_vertical_link_ids = np.empty(shape, np.int)
    node_vertical_link_ids[:-1, :] = vertical_link_ids(shape)
    node_vertical_link_ids[-1, :] = -1

    return node_vertical_link_ids, node_horizontal_link_ids


def _node_in_link_ids(shape):
    """Links entering each node.

    Parameters
    ----------
    shape : tuple of int
        Shape of grid of nodes.

    Returns
    -------
    tuple :
        Tuple of array of link IDs as (vertical_links, horizontal_links).

    Examples
    --------
    >>> from landlab.grid.structured_quad.links import _node_in_link_ids
    >>> (vert, horiz) = _node_in_link_ids((3, 4))
    >>> vert
    array([[-1, -1, -1, -1],
           [ 3,  4,  5,  6],
           [10, 11, 12, 13]])
    >>> horiz
    array([[-1,  0,  1,  2],
           [-1,  7,  8,  9],
           [-1, 14, 15, 16]])
    """
    node_horizontal_link_ids = np.empty(shape, np.int)
    node_horizontal_link_ids[:, 1:] = horizontal_link_ids(shape)
    node_horizontal_link_ids[:, 0] = -1

    node_vertical_link_ids = np.empty(shape, np.int)
    node_vertical_link_ids[1:, :] = vertical_link_ids(shape)
    node_vertical_link_ids[0, :] = -1

    return node_vertical_link_ids, node_horizontal_link_ids


def node_in_link_ids(shape):
    """Links entering each node.

    Parameters
    ----------
    shape : tuple of int
        Shape of grid of nodes.

    Returns
    -------
    tuple :
        Tuple of array of link IDs as (vertical_links, horizontal_links).

    Examples
    --------
    >>> from landlab.grid.structured_quad.links import node_in_link_ids
    >>> (links, offset) = node_in_link_ids((3, 4))
    >>> links
    array([ 0,  1,  2,  3,  4,  7,  5,  8,  6,  9, 10, 11, 14, 12, 15, 13, 16])
    >>> offset
    array([ 0,  0,  1,  2,  3,  4,  6,  8, 10, 11, 13, 15, 17])

    The links entering the 1st, 5th, and last node. The first node does not
    have any links entering it.

    >>> offset[0] == offset[1]
    True
    >>> for link in [4, 11]: links[offset[link]:offset[link + 1]]
    array([3])
    array([13, 16])
    """
    (in_vert, in_horiz) = _node_in_link_ids(shape)
    _node_link_ids = np.vstack((in_vert.flat, in_horiz.flat)).T
    # offset = np.cumsum(number_of_in_links_per_node(shape))

    offset = np.empty(nodes.number_of_nodes(shape) + 1, dtype=int)
    np.cumsum(number_of_in_links_per_node(shape), out=offset[1:])
    offset[0] = 0
    return _node_link_ids[_node_link_ids >= 0], offset


def node_out_link_ids(shape):
    """Links leaving each node.

    Parameters
    ----------
    shape : tuple of int
        Shape of grid of nodes.

    Returns
    -------
    tuple :
        Tuple of array of link IDs as (vertical_links, horizontal_links).

    Examples
    --------
    >>> from landlab.grid.structured_quad.links import node_out_link_ids
    >>> (links, offset) = node_out_link_ids((3, 4))
    >>> links
    array([ 3,  0,  4,  1,  5,  2,  6, 10,  7, 11,  8, 12,  9, 13, 14, 15, 16])
    >>> offset
    array([ 0,  2,  4,  6,  7,  9, 11, 13, 14, 15, 16, 17, 17])

    The links leaving the 1st, 8th, and last node. The last node does not have
    any links leaving it.

    >>> offset[11] == offset[12]
    True
    >>> for link in [0, 7]: links[offset[link]:offset[link + 1]]
    array([3, 0])
    array([13])
    """
    (out_vert, out_horiz) = _node_out_link_ids(shape)
    _node_link_ids = np.vstack((out_vert.flat, out_horiz.flat)).T
    offset = np.empty(nodes.number_of_nodes(shape) + 1, dtype=int)
    np.cumsum(number_of_out_links_per_node(shape), out=offset[1:])
    offset[0] = 0
    return _node_link_ids[_node_link_ids >= 0], offset


def links_at_node(shape):
    """Get link ids for each node.

    Parameters
    ----------
    shape : tuple of int
        Shape of grid of nodes.

    Returns
    -------
    (N, 4) ndarray of int
        Array of link ids.

    Examples
    --------
    >>> from landlab.grid.structured_quad.links import links_at_node
    >>> links_at_node((4, 3)) # doctest: +NORMALIZE_WHITESPACE
    array([[ 0,  2, -1, -1], [ 1,  3,  0, -1], [-1,  4,  1, -1],
           [ 5,  7, -1,  2], [ 6,  8,  5,  3], [-1,  9,  6,  4],
           [10, 12, -1,  7], [11, 13, 10,  8], [-1, 14, 11,  9],
           [15, -1, -1, 12], [16, -1, 15, 13], [-1, -1, 16, 14]])
    """
    (south_links, west_links) = _node_in_link_ids(shape)
    (north_links, east_links) = _node_out_link_ids(shape)

    return np.vstack((east_links.flat, north_links.flat,
                      west_links.flat, south_links.flat)).transpose().copy()


def node_link_ids(shape):
    """Link IDs for links entering and leaving each node.

    Parameters
    ----------
    shape : tuple of int
        Shape of grid of nodes.

    Returns
    -------
    tuple :
        Tuple of array of link IDs and offsets into link array.

    Examples
    --------
    >>> from landlab.grid.structured_quad.links import node_link_ids
    >>> (links, offset) = node_link_ids((3, 4))
    >>> links
    array([ 0,  3,  1,  4,  0,  2,  5,  1,  6,  2,
            7, 10,  3,  8, 11,  7,  4,  9, 12,  8,  5, 13,  9,  6,
           14, 10, 15, 14, 11, 16, 15, 12, 16, 13])
    >>> offset
    array([ 0,  2,  5,  8, 10, 13, 17, 21, 24, 26, 29, 32, 34])

    The links attached to node 0

    >>> links[offset[0]:offset[1]]
    array([0, 3])

    The links attached to node 5

    >>> links[offset[5]:offset[6]]
    array([ 8, 11,  7,  4])
    """
    (in_vert, in_horiz) = _node_in_link_ids(shape)
    (out_vert, out_horiz) = _node_out_link_ids(shape)
    _node_link_ids = np.vstack((out_horiz.flat, out_vert.flat,
                                in_horiz.flat, in_vert.flat)).T

    offset = np.empty(nodes.number_of_nodes(shape) + 1, dtype=int)
    np.cumsum(number_of_links_per_node(shape), out=offset[1:])
    offset[0] = 0

    return _node_link_ids[_node_link_ids >= 0], offset


def node_id_at_link_start(shape):
    """Node ID at start of links.

    Parameters
    ----------
    shape : tuple of int
        Shape of grid of nodes.

    Returns
    -------
    ndarray :
        Node IDs at start of links.

    Examples
    --------
    >>> from landlab.grid.structured_quad.links import node_id_at_link_start
    >>> node_id_at_link_start((3, 4)) # doctest: +NORMALIZE_WHITESPACE
    array([ 0, 1, 2,
            0, 1, 2, 3,
            4, 5, 6,
            4, 5, 6, 7,
            8, 9, 10])
    """
    all_node_ids = nodes.node_ids(shape)
    link_tails_with_extra_row = np.hstack((all_node_ids[:, :-1],
                                           all_node_ids)).reshape((-1, ))
    return link_tails_with_extra_row[:-shape[1]]


def node_id_at_link_end(shape):
    """Node ID at end of links.

    Parameters
    ----------
    shape : tuple of int
        Shape of grid of nodes.

    Returns
    -------
    ndarray :
        Node IDs at end of links.

    Examples
    --------
    >>> from landlab.grid.structured_quad.links import node_id_at_link_end
    >>> node_id_at_link_end((3, 4)) # doctest: +NORMALIZE_WHITESPACE
    array([ 1, 2, 3,
            4, 5, 6, 7,
            5, 6, 7,
            8, 9, 10, 11,
            9, 10, 11])
    """
    all_node_ids = nodes.node_ids(shape)
    link_heads_missing_row = np.hstack((all_node_ids[:-1, 1:],
                                        all_node_ids[1:, :])).reshape((-1, ))
    return np.concatenate((link_heads_missing_row, all_node_ids[-1, 1:]))


def is_active_link(shape, node_status):
    """Link IDs of active links.

    Parameters
    ----------
    shape : tuple of int
        Shape of grid of nodes.
    node_status : array_link
        Status of nodes in grid.

    Returns
    -------
    ndarray :
        Links IDs at the active links.

    Examples
    --------
    >>> from landlab.grid.structured_quad.nodes import (
    ...     status_with_perimeter_as_boundary)
    >>> from landlab.grid.structured_quad.links import is_active_link
    >>> status = status_with_perimeter_as_boundary((3, 4))
    >>> status # doctest: +NORMALIZE_WHITESPACE
    array([[4, 4, 4, 4],
           [4, 0, 0, 4],
           [4, 4, 4, 4]])
    >>> is_active_link((3, 4), status) # doctest: +NORMALIZE_WHITESPACE
    array([False, False, False,
           False, False, False, False,
           False, True, False,
           False, False, False, False,
           False, False, False], dtype=bool)
    """
    if np.prod(shape) != node_status.size:
        raise ValueError('node status array does not match size of grid '
                         '(%d != %d)' % (np.prod(shape), len(node_status)))

    status_at_link_start = node_status.flat[node_id_at_link_start(shape)]
    status_at_link_end = node_status.flat[node_id_at_link_end(shape)]

    return (((status_at_link_start == CORE_NODE) &
             (status_at_link_end == CORE_NODE)) |
            ((status_at_link_end == CORE_NODE) &
             (status_at_link_start == CORE_NODE)) |
            ((status_at_link_end == CORE_NODE) &
             (status_at_link_start == FIXED_VALUE_BOUNDARY)) |
            ((status_at_link_end == FIXED_VALUE_BOUNDARY) &
             (status_at_link_start == CORE_NODE)))


def active_link_ids(shape, node_status):
    """Get active links.

    Parameters
    ----------
    shape : tuple of int
        Shape of grid of nodes.
    node_status : array_link
        Status of nodes in grid.

    Returns
    -------
    ndarray :
        Links IDs at the active links.

    Examples
    --------
    >>> from landlab.grid import RasterModelGrid
    >>> from landlab.grid.structured_quad.links import active_link_ids

    >>> rmg = RasterModelGrid((3, 4))
    >>> rmg.set_closed_boundaries_at_grid_edges(True, True, True, True)

    >>> status = rmg.status_at_node
    >>> status # doctest: +NORMALIZE_WHITESPACE
    array([4, 4, 4, 4,
           4, 0, 0, 4,
           4, 4, 4, 4], dtype=int8)

    >>> active_link_ids((3, 4), status)
    array([8])
    """
    return as_id_array(np.where(is_active_link(shape, node_status))[0])


def is_fixed_link(shape, node_status):
    """ID of active links.

    Parameters
    ----------
    shape : tuple of int
        Shape of grid of nodes.
    node_status : array_link
        Status of nodes in grid.

    Returns
    -------
    ndarray :
        Links IDs at the active links.

    Examples
    --------
    >>> from landlab import RasterModelGrid
    >>> from landlab.grid.structured_quad.links import is_fixed_link
    >>> import numpy as np

    >>> rmg = RasterModelGrid((4, 5))
    >>> z = np.arange(0, rmg.number_of_nodes)
    >>> s = np.arange(0, rmg.number_of_links)
    >>> rmg.at_node['topographic__elevation'] = z
    >>> rmg.at_link['topographic__slope'] = s

    >>> rmg.set_fixed_link_boundaries_at_grid_edges(True, True, True, True)
    >>> rmg.status_at_node # doctest: +NORMALIZE_WHITESPACE
    array([2, 2, 2, 2, 2,
           2, 0, 0, 0, 2,
           2, 0, 0, 0, 2,
           2, 2, 2, 2, 2], dtype=int8)

    >>> is_fixed_link(rmg.shape, rmg.status_at_node)
    array([False, False, False, False, False,  True,  True,  True, False,
            True, False, False,  True, False, False, False, False, False,
            True, False, False,  True, False,  True,  True,  True, False,
           False, False, False, False], dtype=bool)
    """
    if np.prod(shape) != node_status.size:
        raise ValueError('node status array does not match size of grid '
                         '(%d != %d)' % (np.prod(shape), len(node_status)))

    status_at_link_start = node_status.flat[node_id_at_link_start(shape)]
    status_at_link_end = node_status.flat[node_id_at_link_end(shape)]

    return (((status_at_link_start == CORE_NODE) &
             (status_at_link_end == FIXED_GRADIENT_BOUNDARY)) |
            ((status_at_link_end == CORE_NODE) &
             (status_at_link_start == FIXED_GRADIENT_BOUNDARY)))


def fixed_link_ids(shape, node_status):
    """ID of fixed links.

    Parameters
    ----------
    shape : tuple of int
        Shape of grid of nodes.
    node_status : array_link
        Status of nodes in grid.

    Returns
    -------
    ndarray :
        Links IDs at the active links.

    Examples
    --------

    >>> from landlab import RasterModelGrid
    >>> from landlab.grid.structured_quad.links import fixed_link_ids
    >>> import numpy as np
    >>> rmg = RasterModelGrid(4, 5)
    >>> z = np.arange(0, rmg.number_of_nodes)
    >>> s = np.arange(0, rmg.number_of_links)
    >>> rmg.at_node['topographic__elevation'] = z
    >>> rmg.at_link['topographic__slope'] = s
    >>> rmg.set_fixed_link_boundaries_at_grid_edges(True, True, True, True)
    >>> rmg.status_at_node # doctest: +NORMALIZE_WHITESPACE
    array([2, 2, 2, 2, 2,
           2, 0, 0, 0, 2,
           2, 0, 0, 0, 2,
           2, 2, 2, 2, 2], dtype=int8)
    >>> fixed_link_ids(rmg.shape, rmg.status_at_node)
    array([ 5,  6,  7,  9, 12, 18, 21, 23, 24, 25])
    """
    return as_id_array(np.where(is_fixed_link(shape, node_status))[0])


def horizontal_active_link_ids(shape, active_ids, bad_index_value=-1):
    """ID of horizontal active links.

    Parameters
    ----------
    shape : tuple of int
        Shape of grid of nodes.
    active_ids : array of int
        Array of all active link ids
    bad_index_value: int, optional
        Value assigned to inactive indicies in the array.

    Returns
    -------
    ndarray :
        Link IDs at the HORIZONTAL active links. Length of
        number_of_horizontal_links.

    Examples
    --------

    The following example uses this grid::

          *---I-->*---I-->*---I-->*---I-->*
          ^       ^       ^       ^       ^
          I       I       I       I       I
          |       |       |       |       |
          *---I-->o--24-->o--25-->o---I-->*
          ^       ^       ^       ^       ^
          I       V       V       V       I
          |       |       |       |       |
          *---I-->o--20-->o--21-->o---I-->*
          ^       ^       ^       ^       ^
          I       I       I       I       I
          |       |       |       |       |
          *---I-->*---I-->*---I-->*---I-->*

    .. note::

        ``*`` indicates the nodes that are set to :any:`CLOSED_BOUNDARY`

        ``o`` indicates the nodes that are set to :any:`CORE_NODE`

        ``I`` indicates the links that are set to :any:`INACTIVE_LINK`

        ``V`` indicates vertical active ids, which are ignored by this
        function.

        Numeric values correspond to the horizontal :any:`ACTIVE_LINK`  ID.

    >>> from landlab import RasterModelGrid
    >>> from landlab.grid.structured_quad.links import (active_link_ids,
    ...     horizontal_active_link_ids)

    >>> rmg = RasterModelGrid(4, 5)
    >>> rmg.set_closed_boundaries_at_grid_edges(True, True, True, True)

    >>> status = rmg.status_at_node
    >>> status # doctest: +NORMALIZE_WHITESPACE
    array([4, 4, 4, 4, 4,
           4, 0, 0, 0, 4,
           4, 0, 0, 0, 4,
           4, 4, 4, 4, 4], dtype=int8)
    >>> active_ids = active_link_ids((4,5), status)

    >>> horizontal_active_link_ids((4,5), active_ids)
    ...     # doctest: +NORMALIZE_WHITESPACE
    array([-1, -1, -1, -1,
           -1, 10, 11, -1,
           -1, 19, 20, -1,
           -1, -1, -1, -1])
    """
    out = np.full(number_of_horizontal_links(shape), bad_index_value,
                  dtype=int)
    horizontal_ids = active_ids[np.where(~ is_vertical_link(shape, active_ids))]

    out[nth_horizontal_link(shape, horizontal_ids)] = horizontal_ids
    return out


def horizontal_fixed_link_ids(shape, fixed_ids, bad_index_value=-1):
    """ID of horizontal fixed links.

    Parameters
    ----------
    shape : tuple of int
        Shape of grid of nodes.
    fixed_ids : array of int
        Array of all fixed link ids
    bad_index_value: int, optional
        Value assigned to inactive indicies in the array.

    Returns
    -------
    ndarray :
        Link IDs at the HORIZONTAL fixed links. Length of
        number_of_horizontal_links.

    Examples
    --------
    The following example uses this grid::

          *---I-->*---I-->*---I-->*---I-->*
          ^       ^       ^       ^       ^
          I       V       V       V       I
          |       |       |       |       |
          *--18-->o------>o------>o--21-->*
          ^       ^       ^       ^       ^
          I       V       V       V       I
          |       |       |       |       |
          *---9-->o------>o------>o--12-->*
          ^       ^       ^       ^       ^
          I       V       V       V       I
          |       |       |       |       |
          *---I-->*---I-->*---I-->*---I-->*

    .. note::

        ``*`` indicates the nodes that are set to :any:`FIXED_VALUE_BOUNDARY`

        ``o`` indicates the nodes that are set to :any:`CORE_NODE`

        ``I`` indicates the links that are set to :any:`INACTIVE_LINK`

        ``V`` indicates vertical ids, which are ignored by this function

        ``H`` indicates horizontal :any:`ACTIVE_LINK` ids, which are ignored by
        this function

        Numeric values correspond to the horizontal :any:`FIXED_LINK` ID.

    >>> from landlab import RasterModelGrid
    >>> from landlab.grid.structured_quad.links import (fixed_link_ids,
    ...     horizontal_fixed_link_ids)
    >>> import numpy

    >>> rmg = RasterModelGrid(4, 5)
    >>> rmg.at_node['topographic__elevation'] = numpy.arange(
    ...     0, rmg.number_of_nodes)
    >>> rmg.at_link['topographic__slope'] = numpy.arange(
    ...     0, rmg.number_of_links)

    >>> rmg.set_fixed_link_boundaries_at_grid_edges(True, True, True, True)
    >>> status = rmg.status_at_node
    >>> status # doctest: +NORMALIZE_WHITESPACE
    array([2, 2, 2, 2, 2,
           2, 0, 0, 0, 2,
           2, 0, 0, 0, 2,
           2, 2, 2, 2, 2], dtype=int8)

    >>> fixed_ids = fixed_link_ids((4, 5), status)
    >>> horizontal_fixed_link_ids((4, 5), fixed_ids)
    ...     # doctest: +NORMALIZE_WHITESPACE
    array([-1, -1, -1, -1,
            9, -1, -1, 12,
           18, -1, -1, 21,
           -1, -1, -1, -1])
    """
    out = np.full(number_of_horizontal_links(shape), bad_index_value,
                  dtype=int)
    horizontal_ids = fixed_ids[np.where(~ is_vertical_link(shape, fixed_ids))]

    out[nth_horizontal_link(shape, horizontal_ids)] = horizontal_ids
    return out


def is_vertical_link(shape, links):
    """Test if links are vertical.

    Parameters
    ----------
    shape : tuple of int
        Shape of grid of nodes.
    links : array of int
        Array of link ids to test.

    Returns
    -------
    ndarray of bool
        `True` for links that are vertical.

    Examples
    --------
    >>> from landlab.grid.structured_quad.links import (is_vertical_link,
    ...     number_of_links)
    >>> import numpy as np
    >>> shape = (3, 4)
    >>> links = np.arange(number_of_links(shape))
    >>> is_vertical_link(shape, links) # doctest: +NORMALIZE_WHITESPACE
    array([False, False, False,  True,  True,  True,  True,
           False, False, False,  True,  True,  True,  True,
           False, False, False], dtype=bool)
    """
    return (((links % (2 * shape[1] - 1)) >= shape[1] - 1) &
            (links < number_of_links(shape)))


def is_horizontal_link(shape, links):
    """Test if a link is horizontal.

    Parameters
    ----------
    shape : tuple of int
        Shape of grid of nodes.
    links : array of int
        Array of link ids to test.

    Returns
    -------
    ndarray of bool
        `True` for links that are horizontal.

    Examples
    --------
    >>> from landlab.grid.structured_quad.links import (is_horizontal_link,
    ...     number_of_links)
    >>> import numpy as np
    >>> shape = (3, 4)
    >>> links = np.arange(number_of_links(shape))
    >>> is_horizontal_link(shape, links) # doctest: +NORMALIZE_WHITESPACE
    array([ True,  True,  True, False, False, False, False,
            True,  True,  True, False, False, False, False,
            True,  True,  True], dtype=bool)
    """
    return ((~ is_vertical_link(shape, links)) &
            (links < number_of_links(shape)))


def is_diagonal_link(shape, links):
    """Test if a link is diagonal.

    Parameters
    ----------
    shape : tuple of int
        Shape of grid of nodes.
    links : array of int
        Array of link ids to test.

    Returns
    -------
    ndarray of bool
        `True` for links that are diagonal.

    Examples
    --------
    >>> from landlab.grid.structured_quad.links import (is_diagonal_link,
    ...     number_of_links)
    >>> import numpy as np
    >>> shape = (3, 4)
    >>> links = np.array([0, 3, 16, 17])
    >>> is_diagonal_link(shape, links) # doctest: +NORMALIZE_WHITESPACE
    array([False, False, False,  True], dtype=bool)
    """
    return links >= number_of_links(shape)


def nth_vertical_link(shape, links):
    """Convert link ID to vertical link ID.

    Parameters
    ----------
    shape : tuple of int
        Shape of grid of nodes.
    links : array of int
        Array of link ids to test.

    Returns
    -------
    ndarray of int
        The link ID as the nth vertical links.

    Examples
    --------
    >>> from landlab.grid.structured_quad.links import nth_vertical_link
    >>> shape = (3, 4)
    >>> nth_vertical_link(shape, 4)
    1
    >>> nth_vertical_link(shape, (3, 4, 11))
    array([0, 1, 5])
    """
    links = np.asarray(links, dtype=np.int)
    return as_id_array((links // (2 * shape[1] - 1)) * shape[1] +
                       links % (2 * shape[1] - 1) - (shape[1] - 1))


def nth_horizontal_link(shape, links):
    """Convert link ID to horizontal link ID.

    Parameters
    ----------
    shape : tuple of int
        Shape of grid of nodes.
    links : array of int
        Array of link ids to test.

    Returns
    -------
    ndarray of int
        The link ID as the nth horizontal links.

    Examples
    --------
    >>> from landlab.grid.structured_quad.links import nth_horizontal_link
    >>> shape = (3, 4)
    >>> nth_horizontal_link(shape, 16)
    8
    >>> nth_horizontal_link(shape, (1, 7, 8))
    array([1, 3, 4])
    """
    links = np.asarray(links, dtype=np.int)
    return as_id_array((links // (2 * shape[1] - 1)) * (shape[1] - 1) +
                       links % (2 * shape[1] - 1))


def vertical_active_link_ids(shape, active_ids, bad_index_value=-1):
    """ID of vertical active links.

    Parameters
    ----------
    shape : tuple of int
        Shape of grid of nodes.
    active_ids : array of int
        Array of all active link ids
    bad_index_value: int, optional
        Value assigned to inactive indicies in the array.

    Returns
    -------
    ndarray :
        Link IDs at the VERTICAL active links. Length of
        number_of_vertical_links.

    Examples
    --------
    The following example uses this grid::

          *---I-->*---I-->*---I-->*---I-->*
          ^       ^       ^       ^       ^
          I       I       I       I       I
          |       |       |       |       |
          *---I-->o---H-->o---H-->o---I-->*
          ^       ^       ^       ^       ^
          I       6       7       8       I
          |       |       |       |       |
          *---I-->o---H-->o---H-->o---I-->*
          ^       ^       ^       ^       ^
          I       I       I       I       I
          |       |       |       |       |
          *---I-->*---I-->*---I-->*---I-->*

    .. note::

        ``*`` indicates the nodes that are set to :any:`CLOSED_BOUNDARY`

        ``o`` indicates the nodes that are set to :any:`CORE_NODE`

        ``I`` indicates the links that are set to :any:`INACTIVE_LINK`

        ``H`` indicates horizontal active ids, which are ignored by this
        function

        Numeric values correspond to the vertical :any:`ACTIVE_LINK` IDs.

    >>> from landlab import RasterModelGrid
    >>> from landlab.grid.structured_quad.links import (active_link_ids,
    ...     vertical_active_link_ids)

    >>> rmg = RasterModelGrid((4, 5))
    >>> active_ids = active_link_ids((4, 5), rmg.status_at_node)
    >>> active_ids # doctest: +NORMALIZE_WHITESPACE
    array([ 5,  6,  7,
            9, 10, 11, 12,
           14, 15, 16,
           18, 19, 20, 21,
           23, 24, 25])

    >>> vertical_active_link_ids((4, 5), active_ids)
    ...     # doctest: +NORMALIZE_WHITESPACE
    array([-1,  5,  6,  7, -1,
           -1, 14, 15, 16, -1,
           -1, 23, 24, 25, -1])

    >>> rmg.set_closed_boundaries_at_grid_edges(True, True, True, True)
    >>> status = rmg.status_at_node
    >>> active_ids = active_link_ids((4, 5), status)
    >>> vertical_active_link_ids((4, 5), active_ids)
    ...     # doctest: +NORMALIZE_WHITESPACE
    array([-1, -1, -1, -1, -1,
           -1, 14, 15, 16, -1,
           -1, -1, -1, -1, -1])
    """
    out = np.full(number_of_vertical_links(shape), bad_index_value, dtype=int)
    vertical_ids = active_ids[np.where(is_vertical_link(shape, active_ids))]

    out[nth_vertical_link(shape, vertical_ids)] = vertical_ids
    return out


def vertical_fixed_link_ids(shape, fixed_ids, bad_index_value=-1):
    """ID of vertical fixed links.

    Parameters
    ----------
    shape : tuple of int
        Shape of grid of nodes.
    fixed_ids : array of int
        Array of all fixed link ids
    bad_index_value: int, optional
        Value assigned to inactive indicies in the array.

    Returns
    -------
    ndarray :
        Link IDs at the VERTICAL fixed links. Length of
        number_of_vertical_links.

    Examples
    --------
    The following example uses this grid::

        *---I-->*---I-->*---I-->*---I-->*
        ^       ^       ^       ^       ^
        I      23      24      25       I
        |       |       |       |       |
        *---H-->o---H-->o---H-->o---H-->*
        ^       ^       ^       ^       ^
        I       V       V       V       I
        |       |       |       |       |
        *---H-->o---H-->o---H-->o---H-->*
        ^       ^       ^       ^       ^
        I       5       6       7       I
        |       |       |       |       |
        *---I-->*---I-->*---I-->*---I-->*

    .. note::

        ``*`` indicates the nodes that are set to
        :any:`FIXED_GRADIENT_BOUNDARY`

        ``o`` indicates the nodes that are set to :any:`CORE_NODE`

        ``I`` indicates the links that are set to :any:`INACTIVE_LINK`

        ``H`` indicates horizontal active and fixed links, which are ignored by
        this function.

        ``V`` indicates vertical active ids, which are ignored by this
        function.

        Numeric values correspond to the vertical :any:`FIXED_LINK` IDs.

    >>> from landlab import RasterModelGrid
    >>> from landlab.grid.structured_quad.links import (fixed_link_ids,
    ...     vertical_fixed_link_ids)
    >>> import numpy

    >>> rmg = RasterModelGrid((4, 5))
    >>> rmg.at_node['topographic__elevation'] = numpy.arange(
    ...     0, rmg.number_of_nodes)
    >>> rmg.at_link['topographic__slope'] = numpy.arange(
    ...     0, rmg.number_of_links)
    >>> rmg.set_fixed_link_boundaries_at_grid_edges(True, True, True, True)

    >>> status = rmg.status_at_node
    >>> status # doctest: +NORMALIZE_WHITESPACE
    array([2, 2, 2, 2, 2,
           2, 0, 0, 0, 2,
           2, 0, 0, 0, 2,
           2, 2, 2, 2, 2], dtype=int8)
    >>> fixed_ids = fixed_link_ids((4, 5), status)

    >>> vertical_fixed_link_ids((4,5), fixed_ids)
    ...     # doctest: +NORMALIZE_WHITESPACE
    array([-1,  5,  6,  7, -1,
           -1, -1, -1, -1, -1,
           -1, 23, 24, 25, -1])
    """
    out = np.full(number_of_vertical_links(shape), bad_index_value, dtype=int)
    vertical_ids = fixed_ids[np.where(is_vertical_link(shape, fixed_ids))]

    out[nth_vertical_link(shape, vertical_ids)] = vertical_ids
    return out


def horizontal_south_link_neighbor(shape, horizontal_ids,
                                   bad_index_value=-1):
    """ID of south horizontal link neighbor.

    Parameters
    ----------
    shape : tuple of int
        Shape of grid of nodes.
    horizontal_ids : array of int
        Array of all horizontal link ids *must be of len(horizontal_links)*.
    bad_index_value: int, optional
        Value assigned to inactive indicies in the array.

    Returns
    -------
    ndarray :
        Link IDs of south horizontal neighbor links. Length of
        number_of_horizontal_links.

    Examples
    --------
    The following example uses this grid::

        *--27-->*--28-->*--29-->*--30-->*



        *--18-->*--19-->*--20-->*--21-->*



        *---9-->*--10-->*--11-->*--12-->*



        *---0-->*---1-->*---2-->*---3-->*

    .. note::

        Only horizontal links are shown. When no neighbor is found,
        bad_index_value is returned.

        ``*`` indicates nodes

        Numeric values correspond to the horizontal IDs.

    >>> from landlab import RasterModelGrid
    >>> from landlab.grid.structured_quad.links import *
    >>> rmg = RasterModelGrid(4, 5)
    >>> horizontal_links = horizontal_link_ids(rmg.shape).flatten()
    >>> horizontal_south_link_neighbor(rmg.shape, horizontal_links)
    array([-1, -1, -1, -1,  0,  1,  2,  3,  9, 10, 11, 12, 18, 19, 20, 21])
    """
    # First, we find the shape of the horizontal link array given the shape
    # of the raster model grid. In our example, the shape of horizontal links
    # for a grid of 4 rows and 5 columns is 4 rows of horizontal links and 4
    # columns of horizontal links.
    horizontal_2d_shape = shape_of_horizontal_links(shape)

    # Then, we reshape the flattend (1-D) horizontal_link_id array into the
    # shape provided by the shape_of_horizontal_links() function.
    horizontal_2d_array = np.reshape(horizontal_ids, horizontal_2d_shape)

    # To find south links, we need to shift the IDs in the 2-D array. We first
    # insert a row of bad_index_value into the top row of the array
    horizontal_ids = np.insert(horizontal_2d_array, [0], bad_index_value,
                               axis=0)
    # We find the updated array shape and number of rows for the updated array.
    row_len = np.shape(horizontal_ids)[0]

    # To get back to the correct array size (the one found using
    # shape_of_horizontal_links), we delete the last row in the 2-D array
    link_ids = np.delete(horizontal_ids, [row_len - 1], axis=0)

    # Once we have shifted the 2-D array and removed extra indices, we can
    # flatten the output array to a 1-D array with length of
    # number_of_horizontal_links.
    south_horizontal_neighbors = link_ids.flatten()

    return south_horizontal_neighbors


def horizontal_west_link_neighbor(shape, horizontal_ids,
                                  bad_index_value=-1):
    """ID of west, horizontal link neighbor.

    Parameters
    ----------
    shape : tuple of int
        Shape of grid of nodes.
    horizontal_ids : array of int
        Array of all horizontal link ids - *must be of len(horizontal_links)*
    bad_index_value: int, optional
        Value assigned to inactive indicies in the array.

    Returns
    -------
    ndarray :
        Link IDs of west horizontal neighbor links. Length of
        number_of_horizontal_links.

    Examples
    --------
    The following example uses this grid::

        *--27-->*--28-->*--29-->*--30-->*



        *--18-->*--19-->*--20-->*--21-->*



        *---9-->*--10-->*--11-->*--12-->*



        *---0-->*---1-->*---2-->*---3-->*

    .. note::

        Only horizontal links are shown. When no neighbor is found,
        bad_index_value is returned.

        ``*`` indicates nodes

        Numeric values correspond to the horizontal IDs.

    >>> from landlab import RasterModelGrid
    >>> from landlab.grid.structured_quad.links import *
    >>> rmg = RasterModelGrid(4, 5)
    >>> horizontal_links = horizontal_link_ids(rmg.shape).flatten()
    >>> horizontal_west_link_neighbor(rmg.shape, horizontal_links)
    array([-1,  0,  1,  2, -1,  9, 10, 11, -1, 18, 19, 20, -1, 27, 28, 29])
    """
    # First, we find the shape of the horizontal link array given the shape
    # of the raster model grid. In our example, the shape of horizontal links
    # for a grid of 4 rows and 5 columns is 4 rows of horizontal links and 4
    # columns of horizontal links.
    horizontal_2d_shape = shape_of_horizontal_links(shape)

    # Then, we reshape the flattend (1-D) horizontal_link_id array into the
    # shape provided by the shape_of_horizontal_links() function.
    horizontal_2d_array = np.reshape(horizontal_ids, horizontal_2d_shape)

    # To find west links, we need to shift the IDs in the 2-D array. We insert
    # a column of bad_index_value into the first column of the array.
    horizontal_ids = np.insert(horizontal_2d_array, [0], bad_index_value,
                               axis=1)

    # We find the updated array shape and number of columns for the updated
    # array.
    row_len = np.shape(horizontal_ids)[1]

    # To get back to the correct array size (the one found using
    # shape_of_horizontal_links), we delete the very LAST column of the 2-D
    # array. (Any link final column in the 2-D array cannot be a western
    # neighbor anyway).
    horizontal_ids = np.delete(horizontal_ids, [row_len - 1], axis=1)

    # Once we have shifted the 2-D array and removed extra indices, we can
    # flatten the output array to a 1-D array with length of
    # number_of_horizontal_links.
    west_horizontal_neighbors = horizontal_ids.flatten()

    return west_horizontal_neighbors


def horizontal_north_link_neighbor(shape, horizontal_ids,
                                   bad_index_value=-1):
    """ID of north, horizontal link neighbor.

    Parameters
    ----------
    shape : tuple of int
        Shape of grid of nodes.
    horizontal_ids : array of int
        Array of all horizontal link ids - *must be of len(horizontal_links)*
    bad_index_value: int, optional
        Value assigned to inactive indicies in the array.

    Returns
    -------
    ndarray :
        Link IDs of north horizontal neighbor links. Length of
        number_of_horizontal_links.

    Examples
    --------
    The following example uses this grid::

        *--27-->*--28-->*--29-->*--30-->*



        *--18-->*--19-->*--20-->*--21-->*



        *---9-->*--10-->*--11-->*--12-->*



        *---0-->*---1-->*---2-->*---3-->*

    .. note::

        Only horizontal links are shown. When no neighbor is found,
        bad_index_value is returned.

        ``*`` indicates nodes

        Numeric values correspond to the horizontal :any:`ACTIVE_LINK` IDs.

    >>> from landlab import RasterModelGrid
    >>> from landlab.grid.structured_quad.links import *
    >>> rmg = RasterModelGrid(4, 5)
    >>> horizontal_links = horizontal_link_ids(rmg.shape).flatten()
    >>> horizontal_north_link_neighbor(rmg.shape, horizontal_links)
    array([ 9, 10, 11, 12, 18, 19, 20, 21, 27, 28, 29, 30, -1, -1, -1, -1])
    """
    # First, we find the shape of the horizontal link array given the shape
    # of the raster model grid. In our example, the shape of horizontal links
    # for a grid of 4 rows and 5 columns is 4 rows of horizontal links and 4
    # columns of horizontal links.
    horizontal_2d_shape = shape_of_horizontal_links(shape)

    # Then, we reshape the flattend (1-D) horizontal_link_id array into the
    # shape provided by the shape_of_horizontal_links() function.
    horizontal_2d_array = np.reshape(horizontal_ids, horizontal_2d_shape)

    # To find north links, we need to shift the IDs in the 2-D array. We first
    # delete the top row of the array
    horizontal_ids = np.delete(horizontal_2d_array, [0], axis=0)

    # We find the updated array shape and number of rows for the updated array.
    row_len = np.shape(horizontal_ids)[0]

    # To get back to the correct array size (the one found using
    # shape_of_horizontal_links), we insert a row (populated with
    # bad_index_value_ into the end of the 2-D array.
    link_ids = np.insert(horizontal_ids, [row_len], bad_index_value,
                         axis=0)

    # Once we have shifted the 2-D array and removed extra indices, we can
    # flatten the output array to a 1-D array with length of
    # number_of_horizontal_links.
    north_horizontal_neighbors = link_ids.flatten()

    return north_horizontal_neighbors


def horizontal_east_link_neighbor(shape, horizontal_ids,
                                  bad_index_value=-1):
    """IDs of east, horizontal link neighbor.

    Parameters
    ----------
    shape : tuple of int
        Shape of grid of nodes.
    horizontal_ids : array of int
        Array of all horizontal link ids - *must be of len(horizontal_links)*
    bad_index_value: int, optional
        Value assigned to inactive indicies in the array.

    Returns
    -------
    ndarray :
        Link IDs of east horizontal neighbor links. Length of
        number_of_horizontal_links.

    Examples
    --------
    The following example uses this grid::

        *--27-->*--28-->*--29-->*--30-->*



        *--18-->*--19-->*--20-->*--21-->*



        *---9-->*--10-->*--11-->*--12-->*



        *---0-->*---1-->*---2-->*---3-->*

    .. note::

        Only horizontal links are shown. When no neighbor is found,
        bad_index_value is returned.

        ``*`` indicates nodes

        Numeric values correspond to the horizontal :any:`ACTIVE_LINK` IDs.

    >>> from landlab import RasterModelGrid
    >>> from landlab.grid.structured_quad.links import *
    >>> rmg = RasterModelGrid(4, 5)
    >>> horizontal_links = horizontal_link_ids(rmg.shape).flatten()
    >>> horizontal_east_link_neighbor(rmg.shape, horizontal_links)
    array([ 1,  2,  3, -1, 10, 11, 12, -1, 19, 20, 21, -1, 28, 29, 30, -1])
    """
    # First, we find the shape of the horizontal link array given the shape
    # of the raster model grid. In our example, the shape of horizontal links
    # for a grid of 4 rows and 5 columns is 4 rows of horizontal links and 4
    # columns of horizontal links.
    horizontal_2d_shape = shape_of_horizontal_links(shape)

    # Then, we reshape the flattend (1-D) horizontal_link_id array into the
    # shape provided by the shape_of_horizontal_links() function.
    horizontal_2d_array = np.reshape(horizontal_ids, horizontal_2d_shape)

    # To find west links, we need to shift the IDs in the 2-D array. We first
    # delete the first column of the array (these values can never be east
    # neighbors anyway.)
    horizontal_ids = np.delete(horizontal_2d_array, [0], axis=1)

    # We find the updated array shape and number of columns for the updated
    # array.
    row_len = np.shape(horizontal_ids)[1]

    # To get back to the correct array size (the one found using
    # shape_of_horizontal_links), we insert a column of bad_index_value into
    # the last column spot in the 2-D array.
    link_ids = np.insert(horizontal_ids, [row_len], bad_index_value,
                         axis=1)

    # Once we have shifted the 2-D array and removed extra indices, we can
    # flatten the output array to a 1-D array with length of
    # number_of_horizontal_links.
    east_horizontal_neighbors = link_ids.flatten()

    return east_horizontal_neighbors


def d4_horizontal_link_neighbors(shape, horizontal_ids, bad_index_value=-1):
    """IDs of all 4 horizontal link neighbors.

    Parameters
    ----------
    shape : tuple of int
        Shape of grid of nodes.
    horizontal_ids : array of int
        Array of all horizontal link ids - *must be of len(horizontal_links)*
    bad_index_value: int, optional
        Value assigned to inactive indicies in the array.

    Returns
    -------
    ndarray :
        Array of 4 horizontal link neighbors for a given link ID. Returned in
        [E, N, W, S].

    Examples
    --------
    Sample grid, giving neighbors for link ID 10::

        *------>*------>*------>*------>*



        *------>*--19-->*------>*------>*



        *---9-->*--10-->*--11-->*------>*



        *------>*---1-->*------>*------>*

    .. note::

        Only horizontal links are shown. When no neighbor is found,
        bad_index_value is returned.

        ``*`` indicates nodes

        Numeric values correspond to the horizontal IDs.

    >>> from landlab import RasterModelGrid
    >>> from landlab.grid.structured_quad.links import *
    >>> rmg = RasterModelGrid(4, 5)
    >>> horizontal_links = horizontal_link_ids(rmg.shape).flatten()
    >>> d4_horizontal_link_neighbors(rmg.shape, horizontal_links)
    array([[ 1,  9, -1, -1],
           [ 2, 10,  0, -1],
           [ 3, 11,  1, -1],
           [-1, 12,  2, -1],
           [10, 18, -1,  0],
           [11, 19,  9,  1],
           [12, 20, 10,  2],
           [-1, 21, 11,  3],
           [19, 27, -1,  9],
           [20, 28, 18, 10],
           [21, 29, 19, 11],
           [-1, 30, 20, 12],
           [28, -1, -1, 18],
           [29, -1, 27, 19],
           [30, -1, 28, 20],
           [-1, -1, 29, 21]])
    """
    # First we find *south* neighbors...
    south = horizontal_south_link_neighbor(shape, horizontal_ids,
                                           bad_index_value)

    # Then *west* neighbors...
    west = horizontal_west_link_neighbor(shape, horizontal_ids,
                                         bad_index_value)

    # Then *north* neighbors...
    north = horizontal_north_link_neighbor(shape, horizontal_ids,
                                           bad_index_value)

    # Finally, *east* neighbors...
    east = horizontal_east_link_neighbor(shape, horizontal_ids,
                                         bad_index_value)

    # Combine all 4 neighbor arrays into one large array
    # (4 x len_horizontal_links)
    neighbor_array = np.array([east, north, west, south])

    # Transpose the 4 neighbor arrays into a (len_horizontal_links x 4) array.
    neighbor_array = np.transpose(neighbor_array)

    # Output neighbor array. For each input ID, returns [S,W,N,E]
    return neighbor_array


def d4_horizontal_active_link_neighbors(shape, horizontal_ids,
                                        bad_index_value=-1):
    """returns IDs of all 4 horizontal link neighbors.

    Parameters
    ----------
    shape : tuple of int
        Shape of grid of nodes.
    horizontal_ids : array of int
        Array of all horizontal link ids - *must be of len(horizontal_links)*
    bad_index_value: int, optional
        Value assigned to inactive indicies in the array.

    Returns
    -------
    ndarray
        Array of 4 horizontal link neighbors for a given link ID. Returned in
        [E, N, W, S]. Returns array for only ACTIVE horizontal links.

    Examples
    --------
    Sample grid, giving neighbors for link ID 20::


        *------>*------>*------>*------>*



        *------>*--19-->*--20-->*------>*



        *------>*--10-->*--11-->*------>*



        *------>*------>*------>*------>*


    .. note::

        Only horizontal links are shown. When no neighbor is found,
        bad_index_value is returned.

        ``*`` indicates nodes

        Numeric values correspond to the horizontal :any:`ACTIVE_LINK` IDs.

    >>> from landlab import RasterModelGrid
    >>> from landlab.grid.structured_quad.links import *
    >>> rmg = RasterModelGrid((4, 5))
    >>> rmg.set_closed_boundaries_at_grid_edges(True, True, True, True)
    >>> active_ids = active_link_ids(rmg.shape, rmg.status_at_node)
    >>> horizontal_ids = horizontal_active_link_ids(
    ...     rmg.shape, active_ids)
    >>> d4_horizontal_active_link_neighbors(rmg.shape, horizontal_ids)
    array([[11, 19, -1, -1],
           [-1, 20, 10, -1],
           [20, -1, -1, 10],
           [-1, -1, 19, 11]])
    """
    # To do this we simply call the find_d4_horizontal_neighbors() function
    # which gives the neighbors for ALL horizontal links in an array, even
    # inactive links.
    d4_neigh = d4_horizontal_link_neighbors(shape, horizontal_ids,
                                            bad_index_value)

    # Now we will just focus on indices that are ACTIVE...
    active_links = np.where(horizontal_ids != bad_index_value)

    # Clip our initial array into a smaller one with just active neighbors
    neighbor_array = d4_neigh[active_links]

    # Output neighbor array. For each input ID, returns [S,W,N,E]
    return neighbor_array


def vertical_south_link_neighbor(shape, vertical_ids, bad_index_value=-1):
    """Link IDs of south, vertical link neighbor.

    Parameters
    ----------
    shape : tuple of int
        Shape of grid of nodes.
    vertical_ids : array of int
        Array of all vertical link ids - MUST BE ARRAY OF LEN(VERTICAL_LINKS)
    bad_index_value: int, optional
        Value assigned to inactive indicies in the array.

    Returns
    -------
    ndarray :
        Link IDs of *south* vertical neighbor links. Length of
        number_of_vertical_links.

    Examples
    --------
    The following example uses this grid::

        *       *       *       *       *
        ^       ^       ^       ^       ^
       22      23      24      25      26
        |       |       |       |       |
        *       *       *       *       *
        ^       ^       ^       ^       ^
       13      14      15      16      17
        |       |       |       |       |
        *       *       *       *       *
        ^       ^       ^       ^       ^
        4       5       6       7       8
        |       |       |       |       |
        *       *       *       *       *

    .. note::

        Only vertical links are shown. When no neighbor is found,
        bad_index_value is returned.

        ``*`` indicates nodes

        Numeric values correspond to the vertical IDs.

    >>> from landlab import RasterModelGrid
    >>> from landlab.grid.structured_quad.links import *
    >>> rmg = RasterModelGrid((4, 5))
    >>> vertical_links = vertical_link_ids(rmg.shape)
    >>> vertical_south_link_neighbor(rmg.shape, vertical_links)
    ...     # doctest: +NORMALIZE_WHITESPACE
    array([-1, -1, -1, -1, -1,
            4,  5,  6,  7,  8,
            13, 14, 15, 16, 17])
    """
    # First, we find the shape of the vertical link array given the shape
    # of the raster model grid. In our example, the shape of vertical links for
    # a grid of 4 rows and 5 columns is 3 rows of vertical links and 5 columns
    # of vertical links.
    vertical_2d_shape = shape_of_vertical_links(shape)

    # Then, we reshape the flattend (1-D) vertical_link_id array into the shape
    # provided by the shape_of_vertical_links() function.
    vertical_2d_array = np.reshape(vertical_ids, vertical_2d_shape)

    # To find south links, we need to shift the IDs in the 2-D array. We insert
    # a row of bad_index_value into the top row of the 2-D array
    link_ids = np.insert(vertical_2d_array, [0], bad_index_value, axis=0)

    # We find the updated array shape and number of rows for the updated array.
    row_len = np.shape(link_ids)[0]

    # To get back to the correct array size (the one found using
    # shape_of_vertical_links), we delete a the last row of the 2-D array.
    vertical_ids = np.delete(link_ids, [row_len - 1], axis=0)

    # Once we have shifted the 2-D array and removed extra indices, we can
    # flatten the output array to a 1-D array with length of
    # number_of_vertical_links.
    south_vertical_neighbors = vertical_ids.flatten()

    return south_vertical_neighbors


def vertical_west_link_neighbor(shape, vertical_ids, bad_index_value=-1):
    """Link IDs of west, vertical link neighbor.

    Parameters
    ----------
    shape : tuple of int
        Shape of grid of nodes.
    vertical_ids : array of int
        Array of all vertical link ids- MUST BE ARRAY OF LEN(VERTICAL_LINKS)
    bad_index_value: int, optional
        Value assigned to inactive indicies in the array.

    Returns
    -------
    ndarray :
        Link IDs of *west* vertical neighbor links. Length of
        number_of_vertical_links.

    Examples
    --------
    The following example uses this grid::

        *       *       *       *       *
        ^       ^       ^       ^       ^
       22      23      24      25      26
        |       |       |       |       |
        *       *       *       *       *
        ^       ^       ^       ^       ^
       13      14      15      16      17
        |       |       |       |       |
        *       *       *       *       *
        ^       ^       ^       ^       ^
        4       5       6       7       8
        |       |       |       |       |
        *       *       *       *       *

    .. note::

        Only vertical links are shown. When no neighbor is found,
        bad_index_value is returned.

        ``*`` indicates nodes

        Numeric values correspond to the vertical IDs.

    >>> from landlab import RasterModelGrid
    >>> from landlab.grid.structured_quad.links import *
    >>> rmg = RasterModelGrid(4, 5)
    >>> vertical_links = vertical_link_ids(rmg.shape)
    >>> vertical_west_link_neighbor(rmg.shape, vertical_links)
    ...     # doctest: +NORMALIZE_WHITESPACE
    array([-1,  4,  5,  6,  7,
           -1, 13, 14, 15, 16,
           -1, 22, 23, 24, 25])
    """
    # First, we find the shape of the vertical link array given the shape
    # of the raster model grid. In our example, the shape of vertical links for
    # a grid of 4 rows and 5 columns is 3 rows of vertical links and 5 columns
    # of vertical links.
    vertical_2d_shape = shape_of_vertical_links(shape)

    # Then, we reshape the flattend (1-D) vertical_link_id array into the shape
    # provided by the shape_of_vertical_links() function.
    vertical_2d_array = np.reshape(vertical_ids, vertical_2d_shape)

    # To find west links, we need to shift the IDs in the 2-D array. We insert
    # a column of bad_index_value into the first column of the array.
    vertical_ids = np.insert(vertical_2d_array, [0], bad_index_value,
                             axis=1)

    # We find the updated array shape and number of columns for the updated
    # array.
    row_len = np.shape(vertical_ids)[1]

    # To get back to the correct array size (the one found using
    # shape_of_vertical_links), we delete the very LAST column of the 2-D
    # array. (Any link final column in the 2-D array cannot be a western
    # neighbor anyway).
    vertical_ids = np.delete(vertical_ids, [row_len - 1], axis=1)

    # Once we have shifted the 2-D array and removed extra indices, we can
    # flatten the output array to a 1-D array with length of
    # number_of_vertical_links.
    west_vertical_neighbors = vertical_ids.flatten()

    return west_vertical_neighbors


def vertical_north_link_neighbor(shape, vertical_ids, bad_index_value=-1):
    """Link IDs of north, vertical link neighbor.

    Parameters
    ----------
    shape : tuple of int
        Shape of grid of nodes.
    vertical_ids : array of int
        Array of all vertical link ids- MUST BE ARRAY OF LEN(VERTICAL_LINKS)
    bad_index_value: int, optional
        Value assigned to inactive indicies in the array.

    Returns
    -------
    ndarray :
        Link IDs of *north* vertical neighbor links. Length of
        number_of_vertical_links.

    Examples
    --------
    The following example uses this grid::

        *       *       *       *       *
        ^       ^       ^       ^       ^
       22      23      24      25      26
        |       |       |       |       |
        *       *       *       *       *
        ^       ^       ^       ^       ^
       13      14      15      16      17
        |       |       |       |       |
        *       *       *       *       *
        ^       ^       ^       ^       ^
        4       5       6       7       8
        |       |       |       |       |
        *       *       *       *       *

    .. note::

        Only vertical links are shown. When no neighbor is found,
        bad_index_value is returned.

        ``*`` indicates nodes

        Numeric values correspond to the vertical IDs.

    >>> from landlab import RasterModelGrid
    >>> from landlab.grid.structured_quad.links import *
    >>> rmg = RasterModelGrid(4, 5)
    >>> vertical_ids = vertical_link_ids(rmg.shape)
    >>> vertical_north_link_neighbor(rmg.shape, vertical_ids)
    ...     # doctest: +NORMALIZE_WHITESPACE
    array([13, 14, 15, 16, 17,
           22, 23, 24, 25, 26,
           -1, -1, -1, -1, -1])
    """
    # First, we find the shape of the vertical link array given the shape
    # of the raster model grid. In our example, the shape of vertical links for
    # a grid of 4 rows and 5 columns is 3 rows of vertical links and 5 columns
    # of vertical links.
    vertical_2d_shape = shape_of_vertical_links(shape)

    # Then, we reshape the flattend (1-D) vertical_link_id array into the shape
    # provided by the shape_of_vertical_links() function.
    vertical_2d_array = np.reshape(vertical_ids, vertical_2d_shape)

    # To find north links, we need to shift the IDs in the 2-D array. We first
    # delete the first row of the array.
    vertical_ids = np.delete(vertical_2d_array, [0], axis=0)

    # We find the updated array shape and number of rows for the updated array.
    row_len = np.shape(vertical_ids)[0]

    # To get back to the correct array size (the one found using
    # shape_of_vertical_links), we insert a row (populated with
    # bad_index_value) into the end of the 2-D array.
    link_ids = np.insert(vertical_ids, [row_len], bad_index_value,
                         axis=0)

    # Once we have shifted the 2-D array and removed extra indices, we can
    # flatten the output array to a 1-D array with length of
    # number_of_vertical_links.
    north_vertical_neighbors = link_ids.flatten()

    return north_vertical_neighbors


def vertical_east_link_neighbor(shape, vertical_ids, bad_index_value=-1):
    """Link IDs of east, vertical link neighbor.

    Parameters
    ----------
    shape : tuple of int
        Shape of grid of nodes.
    vertical_ids : array of int
        Array of all vertical link ids - MUST BE ARRAY OF LEN(VERTICAL_LINKS)
    bad_index_value: int, optional
        Value assigned to inactive indicies in the array.

    Returns
    -------
    ndarray :
        Link IDs of *east* vertical neighbor links. Length of
        number_of_vertical_links.

    Examples
    --------
    The following example uses this grid::

        *       *       *       *       *
        ^       ^       ^       ^       ^
       22      23      24      25      26
        |       |       |       |       |
        *       *       *       *       *
        ^       ^       ^       ^       ^
       13      14      15      16      17
        |       |       |       |       |
        *       *       *       *       *
        ^       ^       ^       ^       ^
        4       5       6       7       8
        |       |       |       |       |
        *       *       *       *       *

    .. note::

        Only vertical links are shown. When no neighbor is found,
        bad_index_value is returned.

        ``*`` indicates nodes

        Numeric values correspond to the vertical IDs.

    >>> from landlab import RasterModelGrid
    >>> from landlab.grid.structured_quad.links import *
    >>> rmg = RasterModelGrid(4, 5)
    >>> vertical_links = vertical_link_ids(rmg.shape)
    >>> vertical_east_link_neighbor(rmg.shape, vertical_links)
    ...     # doctest: +NORMALIZE_WHITESPACE
    array([ 5,  6,  7,  8, -1,
           14, 15, 16, 17, -1,
           23, 24, 25, 26, -1])
    """
    # First, we find the shape of the vertical link array given the shape
    # of the raster model grid. In our example, the shape of vertical links for
    # a grid of 4 rows and 5 columns is 3 rows of vertical links and 5 columns
    # of vertical links.
    vertical_2d_shape = shape_of_vertical_links(shape)

    # Then, we reshape the flattend (1-D) vertical_link_id array into the shape
    # provided by the shape_of_vertical_links() function.
    vertical_2d_array = np.reshape(vertical_ids, vertical_2d_shape)

    # To find east links, we need to shift the IDs in the 2-D array. We first
    # delete the first column of the array.
    vertical_ids = np.delete(vertical_2d_array, [0], axis=1)

    # We find the updated array shape and number of columns for the updated
    # array.
    row_len = np.shape(vertical_ids)[1]

    # To get back to the correct array size (the one found using
    # shape_of_vertical_links), we insert a column (populated with
    # bad_index_value) into the end of the 2-D array.
    link_ids = np.insert(vertical_ids, [row_len], bad_index_value,
                         axis=1)

    # Once we have shifted the 2-D array and removed extra indices, we can
    # flatten the output array to a 1-D array with length of
    # number_of_vertical_links.
    east_vertical_neighbors = link_ids.flatten()

    return east_vertical_neighbors


def d4_vertical_link_neighbors(shape, vertical_ids, bad_index_value=-1):
    """IDs of all 4 vertical link neighbors.

    Parameters
    ----------
    shape : tuple of int
        Shape of grid of nodes.
    vertical_ids : array of int
        Array of all vertical link ids - MUST BE ARRAY OF LEN(VERTICAL_LINKS)
    bad_index_value: int, optional
        Value assigned to inactive indicies in the array.

    Returns
    -------
    ndarray :
        Array of 4 vertical link neighbors for a given link ID. Returned in
        [E, N, W, S].

    Examples
    --------
    The following example uses this grid::

        *       *       *       *       *
        ^       ^       ^       ^       ^
       22      23      24      25      26
        |       |       |       |       |
        *       *       *       *       *
        ^       ^       ^       ^       ^
       13      14      15      16      17
        |       |       |       |       |
        *       *       *       *       *
        ^       ^       ^       ^       ^
        4       5       6       7       8
        |       |       |       |       |
        *       *       *       *       *

    .. note::

        Only vertical links are shown. When no neighbor is found,
        bad_index_value is returned.

        ``*`` indicates nodes

        Numeric values correspond to the vertical IDs.

    >>> from landlab import RasterModelGrid
    >>> from landlab.grid.structured_quad.links import *
    >>> rmg = RasterModelGrid(4, 5)
    >>> vertical_ids = vertical_link_ids(rmg.shape)
    >>> d4_vertical_link_neighbors(rmg.shape, vertical_ids)
    array([[ 5, 13, -1, -1],
           [ 6, 14,  4, -1],
           [ 7, 15,  5, -1],
           [ 8, 16,  6, -1],
           [-1, 17,  7, -1],
           [14, 22, -1,  4],
           [15, 23, 13,  5],
           [16, 24, 14,  6],
           [17, 25, 15,  7],
           [-1, 26, 16,  8],
           [23, -1, -1, 13],
           [24, -1, 22, 14],
           [25, -1, 23, 15],
           [26, -1, 24, 16],
           [-1, -1, 25, 17]])
    """
    south = vertical_south_link_neighbor(shape, vertical_ids, bad_index_value)
    west = vertical_west_link_neighbor(shape, vertical_ids, bad_index_value)
    north = vertical_north_link_neighbor(shape, vertical_ids, bad_index_value)
    east = vertical_east_link_neighbor(shape, vertical_ids, bad_index_value)
    neighbor_array = np.array([east, north, west, south])
    neighbor_array = np.transpose(neighbor_array)
    return neighbor_array


def d4_vertical_active_link_neighbors(shape, vertical_ids, bad_index_value=-1):
    """IDs of all 4 vertical link neighbors.

    Parameters
    ----------
    shape : tuple of int
        Shape of grid of nodes.
    vertical_active_ids : array of int
        Array of all vertical link ids - *must be of len(vertical_links)*
    bad_index_value: int, optional
        Value assigned to inactive indicies in the array.

    Returns
    -------
    ndarray :
        Array of 4 vertical link neighbors for a given ACTIVE link ID.
        Returned in [E, N, W, S].

    Examples
    --------
    The following example uses this grid::

        *       *       *       *       *
        ^       ^       ^       ^       ^
       22      23      24      25      26
        |       |       |       |       |
        *       *       *       *       *
        ^       ^       ^       ^       ^
       13      14      15      16      17
        |       |       |       |       |
        *       *       *       *       *
        ^       ^       ^       ^       ^
        4       5       6       7       8
        |       |       |       |       |
        *       *       *       *       *

    .. note::

        Only vertical links are shown. When no neighbor is found,
        bad_index_value is returned.

        ``*`` indicates nodes

        Numeric values correspond to the vertical IDs.

    >>> from landlab import RasterModelGrid
    >>> from landlab.grid.structured_quad.links import (active_link_ids,
    ...     vertical_active_link_ids, d4_vertical_active_link_neighbors)
    >>> rmg = RasterModelGrid(4, 5)
    >>> active_link_ids = active_link_ids(rmg.shape, rmg.status_at_node)
    >>> vertical_active_ids = vertical_active_link_ids(
    ...     rmg.shape, active_link_ids)
    >>> d4_vertical_active_link_neighbors(rmg.shape, vertical_active_ids)
    array([[ 6, 14, -1, -1],
           [ 7, 15,  5, -1],
           [-1, 16,  6, -1],
           [15, 23, -1,  5],
           [16, 24, 14,  6],
           [-1, 25, 15,  7],
           [24, -1, -1, 14],
           [25, -1, 23, 15],
           [-1, -1, 24, 16]])
    """
    # To do this we simply call the find_d4_vertical_neighbors() function
    # which gives the neighbors for ALL vertical links in an array, even
    # inactive links.
    d4_all_neighbors = d4_vertical_link_neighbors(shape, vertical_ids,
                                                  bad_index_value)

    # Now we will just focus on indices that are ACTIVE...
    active_links = np.where(vertical_ids != bad_index_value)

    # Clip our initial array into a smaller one with just active neighbors
    neighbor_array = d4_all_neighbors[active_links]

    # Output neighbor array. For each input ID, returns [S,W,N,E]
    return neighbor_array


def bottom_edge_horizontal_ids(shape):
    """Link IDs of bottom edge horizontal links.

    Parameters
    ----------
    shape : tuple of int
        Shape of grid, given as (rows, columns) of nodes.

    Returns
    -------
    ndarray :
        Link IDs of bottom edge horizontal links. Length is
        (rmg.number_of_columns-1)

    Examples
    --------
    The following example uses this grid::

        *--27-->*--28-->*--29-->*--30-->*



        *--18-->*--19-->*--20-->*--21-->*



        *---9-->*--10-->*--11-->*--12-->*



        *---0-->*---1-->*---2-->*---3-->*

    .. note::

        Only horizontal links are shown.

        ``*`` indicates nodes

        Numeric values correspond to the horizontal IDs.

    >>> from landlab import RasterModelGrid
    >>> from landlab.grid.structured_quad.links import (
    ...     bottom_edge_horizontal_ids)
    >>> rmg = RasterModelGrid(4, 5)
    >>> shape = rmg.shape
    >>> bottom_edge_horizontal_ids(shape)
    array([0, 1, 2, 3])
    """
    # First, we find all horizontal link ids for the RasterModelGrid shape.
    horizontal_id_array = horizontal_link_ids(shape)

    # Then we slice the first column and return it. This has our bottom edge
    # horizontal ids. This array should be equal in length to (number of
    # columns - 1)
    bottom_edge_hori_ids = horizontal_id_array[0]

    return bottom_edge_hori_ids


def left_edge_horizontal_ids(shape):
    """Link IDs of left edge horizontal links.

    Parameters
    ----------
    shape : tuple of int
        Shape of grid, given as (rows, columns) of nodes.

    Returns
    -------
    ndarray :
        Link IDs of left edge horizontal links. Length is (rmg.number_of_rows)

    Examples
    --------
    The following example uses this grid::

        *--27-->*--28-->*--29-->*--30-->*



        *--18-->*--19-->*--20-->*--21-->*



        *---9-->*--10-->*--11-->*--12-->*



        *---0-->*---1-->*---2-->*---3-->*

    .. note::

        Only horizontal links are shown.

        ``*`` indicates nodes

        Numeric values correspond to the horizontal IDs.

    >>> from landlab import RasterModelGrid
    >>> from landlab.grid.structured_quad.links import left_edge_horizontal_ids
    >>> rmg = RasterModelGrid(4, 5)
    >>> shape = rmg.shape
    >>> left_edge_horizontal_ids(shape)
    array([ 0,  9, 18, 27])
    """

    # First, we find all horizontal link ids for the RasterModelGrid shape.
    horizontal_id_array = horizontal_link_ids(shape)

    # Then we slice the first column and return it. This has our left edge
    # horizontal ids. This array should be equal in length to (number of rows)
    left_edge_hori_ids = horizontal_id_array[:, 0]

    return left_edge_hori_ids


def top_edge_horizontal_ids(shape):
    """IDs of top edge horizontal links.

    Parameters
    ----------
    shape : tuple of int
        Shape of grid, given as (rows, columns) of nodes.

    Returns
    -------
    ndarray :
        Link IDs of top edge horizontal links. Length is
        (rmg.number_of_columns - 1)

    Examples
    --------
    The following example uses this grid::

        *--27-->*--28-->*--29-->*--30-->*



        *--18-->*--19-->*--20-->*--21-->*



        *---9-->*--10-->*--11-->*--12-->*



        *---0-->*---1-->*---2-->*---3-->*

    .. note::

        Only horizontal links are shown.

        ``*`` indicates nodes

        Numeric values correspond to the horizontal IDs.

    >>> from landlab import RasterModelGrid
    >>> from landlab.grid.structured_quad.links import top_edge_horizontal_ids
    >>> rmg = RasterModelGrid(4, 5)
    >>> shape = rmg.shape
    >>> top_edge_horizontal_ids(shape)
    array([27, 28, 29, 30])
    """
    # First, we find all horizontal link ids for the RasterModelGrid shape.
    horizontal_id_array = horizontal_link_ids(shape)

    # Then we slice the first column and return it. This has our top edge
    # horizontal ids. This array should be equal in length to (number of
    # columns - 1)
    top_edge_hori_ids = horizontal_id_array[(shape[0] - 1)]

    return top_edge_hori_ids


def right_edge_horizontal_ids(shape):
    """IDs of right edge horizontal links.

    Parameters
    ----------
    shape : tuple of int
        Shape of grid, given as (rows, columns) of nodes.

    Returns
    -------
    ndarray :
        Link IDs of left edge horizontal links. Length is (rmg.number_of_rows)

    Examples
    --------
    The following example uses this grid::

        *--27-->*--28-->*--29-->*--30-->*



        *--18-->*--19-->*--20-->*--21-->*



        *---9-->*--10-->*--11-->*--12-->*



        *---0-->*---1-->*---2-->*---3-->*

    .. note::

        Only horizontal links are shown.

        ``*`` indicates nodes

        Numeric values correspond to the horizontal IDs.

    >>> from landlab import RasterModelGrid
    >>> from landlab.grid.structured_quad.links import (
    ...     right_edge_horizontal_ids)
    >>> rmg = RasterModelGrid(4, 5)
    >>> shape = rmg.shape
    >>> right_edge_horizontal_ids(shape)
    array([ 3, 12, 21, 30])
    """
    # First, we find all horizontal link ids for the RasterModelGrid shape.
    horizontal_id_array = horizontal_link_ids(shape)

    # Then we slice the last column and return it. This has our right edge
    # horizontal ids. This array should be equal in length to (number of
    # columns - 2)
    right_edge_hori_ids = horizontal_id_array[:, (shape[1] - 2)]

    return right_edge_hori_ids


def bottom_edge_vertical_ids(shape):
    """Link IDs of bottom edge vertical links.

    Parameters
    ----------
    shape : tuple of int
        Shape of grid, given as (rows, columns) of nodes.

    Returns
    -------
    ndarray :
        Link IDs of bottom edge vertical links. Length is
        (rmg.number_of_columns)

    Examples
    --------
    The following example uses this grid::

        *       *       *       *       *
        ^       ^       ^       ^       ^
       22      23      24      25      26
        |       |       |       |       |
        *       *       *       *       *
        ^       ^       ^       ^       ^
       13      14      15      16      17
        |       |       |       |       |
        *       *       *       *       *
        ^       ^       ^       ^       ^
        4       5       6       7       8
        |       |       |       |       |
        *       *       *       *       *

    .. note::

        Only vertical links are shown.

        ``*`` indicates nodes

        Numeric values correspond to the vertical IDs.

    >>> from landlab import RasterModelGrid
    >>> from landlab.grid.structured_quad.links import bottom_edge_vertical_ids
    >>> rmg = RasterModelGrid(4, 5)
    >>> shape = rmg.shape
    >>> bottom_edge_vertical_ids(shape)
    array([4, 5, 6, 7, 8])
    """
    # First, we find all vertical link ids for the RasterModelGrid shape.
    vertical_id_array = vertical_link_ids(shape)

    # Then we slice the first column and return it. This has our bottom edge
    # vertical ids. This array should be equal in length to (number of columns)
    bottom_edge_vert_ids = vertical_id_array[0]

    return bottom_edge_vert_ids


def left_edge_vertical_ids(shape):
    """Link IDs of left edge vertical links.

    Parameters
    ----------
    shape : tuple of int
        Shape of grid, given as (rows, columns) of nodes.

    Returns
    -------
    ndarray :
        Link IDs of left edge vertical links. Length is
        (rmg.number_of_rows - 1)

    Examples
    --------
    The following example uses this grid::

        *       *       *       *       *
        ^       ^       ^       ^       ^
       22      23      24      25      26
        |       |       |       |       |
        *       *       *       *       *
        ^       ^       ^       ^       ^
       13      14      15      16      17
        |       |       |       |       |
        *       *       *       *       *
        ^       ^       ^       ^       ^
        4       5       6       7       8
        |       |       |       |       |
        *       *       *       *       *

    .. note::

        Only vertical links are shown.

        ``*`` indicates nodes

        Numeric values correspond to the vertical IDs.

    >>> from landlab import RasterModelGrid
    >>> from landlab.grid.structured_quad.links import left_edge_vertical_ids
    >>> rmg = RasterModelGrid(4, 5)
    >>> shape = rmg.shape
    >>> left_edge_vertical_ids(shape)
    array([ 4, 13, 22])
    """
    # First, we find all vertical link ids for the RasterModelGrid shape.
    vertical_id_array = vertical_link_ids(shape)

    # Then we slice the first column and return it. This has our left edge
    # vertical ids. This array should be equal in length to
    # (number of rows - 1)
    left_edge_vert_ids = vertical_id_array[:, 0]

    return left_edge_vert_ids


def top_edge_vertical_ids(shape):
    """Link IDs of top edge vertical links.

    Parameters
    ----------
    shape : tuple of int
        Shape of grid, given as (rows, columns) of nodes.

    Returns
    -------
    ndarray :
        Link IDs of top edge vertical links. Length is (rmg.number_of_columns)

    Examples
    --------
    The following example uses this grid::

        *       *       *       *       *
        ^       ^       ^       ^       ^
       22      23      24      25      26
        |       |       |       |       |
        *       *       *       *       *
        ^       ^       ^       ^       ^
       13      14      15      16      17
        |       |       |       |       |
        *       *       *       *       *
        ^       ^       ^       ^       ^
        4       5       6       7       8
        |       |       |       |       |
        *       *       *       *       *

    .. note::

        Only vertical links are shown.

        ``*`` indicates nodes

        Numeric values correspond to the vertical IDs.

    >>> from landlab import RasterModelGrid
    >>> from landlab.grid.structured_quad.links import top_edge_vertical_ids
    >>> rmg = RasterModelGrid(4, 5)
    >>> shape = rmg.shape
    >>> top_edge_vertical_ids(shape)
    array([22, 23, 24, 25, 26])
    """
    # First, we find all vertical link ids for the RasterModelGrid shape.
    vertical_id_array = vertical_link_ids(shape)

    # Then we slice the first column and return it. This has our top edge
    # vertical ids. This array should be equal in length to (number of columns)
    top_edge_vert_ids = vertical_id_array[(shape[0] - 2)]

    return top_edge_vert_ids


def right_edge_vertical_ids(shape):
    """Link IDs of right edge vertical links.

    Parameters
    ----------
    shape : tuple of int
        Shape of grid, given as (rows, columns) of nodes.

    Returns
    -------
    ndarray :
        Link IDs of left edge vertical links. Length is
        (rmg.number_of_rows - 1)

    Examples
    --------
    The following example uses this grid::

        *       *       *       *       *
        ^       ^       ^       ^       ^
       22      23      24      25      26
        |       |       |       |       |
        *       *       *       *       *
        ^       ^       ^       ^       ^
       13      14      15      16      17
        |       |       |       |       |
        *       *       *       *       *
        ^       ^       ^       ^       ^
        4       5       6       7       8
        |       |       |       |       |
        *       *       *       *       *

    .. note::

        Only vertical links are shown.

        ``*`` indicates nodes

        Numeric values correspond to the vertical IDs.

    >>> from landlab import RasterModelGrid
    >>> from landlab.grid.structured_quad.links import right_edge_vertical_ids
    >>> rmg = RasterModelGrid(4, 5)
    >>> shape = rmg.shape
    >>> right_edge_vertical_ids(shape)
    array([ 8, 17, 26])
    """
    # First, we find all vertical link ids for the RasterModelGrid shape.
    vertical_id_array = vertical_link_ids(shape)

    # Then we slice the last column and return it. This has our right edge
    # vertical ids. This array should be equal in length to
    # (number of rows - 1)
    right_edge_vert_ids = vertical_id_array[:, (shape[1] - 1)]

    return right_edge_vert_ids


class StructuredQuadLinkGrid(LinkGrid):

    def __init__(self, shape):
        link_ends = (node_id_at_link_start(shape), node_id_at_link_end(shape))
        number_of_nodes = np.prod(shape)
        LinkGrid.__init__(self, link_ends, number_of_nodes)

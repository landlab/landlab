#! /usr/bin/env python
"""Utility functions for structured grid of elements with four neighbors."""


import itertools

import numpy as np
from six.moves import range

from ..core.utils import as_id_array
from ..grid.base import (
    BAD_INDEX_VALUE,
    CLOSED_BOUNDARY,
    CORE_NODE,
    FIXED_VALUE_BOUNDARY,
)


def node_count(shape):
    """Total number of nodes.

    The total number of nodes in a structured grid with dimensions given
    by the tuple, *shape*. Where *shape* is the number of node rows and
    node columns.

    >>> from landlab.utils.structured_grid import node_count
    >>> node_count((3, 4))
    12
    """
    assert len(shape) == 2
    return shape[0] * shape[1]


def interior_node_count(shape):
    """Number of interior nodes.

    Return the count of the number of interior nodes of a structured grid
    of dimensions, *shape*.

    >>> from landlab.utils.structured_grid import (node_count,
    ...                                            interior_node_count)
    >>> node_count((2, 4))
    8
    >>> interior_node_count((2, 4))
    0
    >>> interior_node_count((1, 4))
    0
    >>> interior_node_count((3, 4))
    2
    """
    assert len(shape) == 2
    if np.min(shape) > 2:
        return (shape[0] - 2) * (shape[1] - 2)
    else:
        return 0


def cell_count(shape):
    """Total number of cells.

    The total number of cells in a structured grid with dimensions, *shape*.
    Where *shape* is a tuple that gives the dimensions of the grid as number
    of rows of nodes followed by number of columns of nodes.

    >>> from landlab.utils.structured_grid import cell_count
    >>> cell_count((3, 4))
    2
    >>> cell_count((1, 4))
    0
    """
    assert len(shape) == 2
    if np.min(shape) > 2:
        return (shape[0] - 2) * (shape[1] - 2)
    else:
        return 0


def active_cell_count(shape):
    """Number of active cells.

    Number of active cells. By default, all cells are active so this is
    the same as cell_count. (active = core+open boundary)
    """
    return cell_count(shape)


def core_cell_count(shape):
    """Number of core cells.

    Number of core cells. By default, all cells are core so this is
    the same as cell_count.
    """
    return cell_count(shape)


def active_link_count(shape):
    """Number of active links.

    Number of active links in a structured grid with dimensions, *shape*.
    A link is active if it connects to at least one active node.

    >>> from landlab.utils.structured_grid import link_count, active_link_count
    >>> link_count((3, 2))
    7
    >>> active_link_count((3, 2))
    0
    >>> active_link_count((3, 4))
    7
    """
    assert len(shape) == 2
    if np.min(shape) > 2:
        return 2 * shape[0] * shape[1] - 3 * (shape[0] + shape[1]) + 4
    else:
        return 0


def link_count(shape):
    """Total number of links.

    Total (active and inactive) number of links in a structured grid with
    dimensions, *shape*. This is the number of to-links and from-links, not
    the total of the two.

    >>> from landlab.utils.structured_grid import link_count
    >>> link_count((3,2))
    7
    """
    assert len(shape) == 2
    return shape[1] * (shape[0] - 1) + shape[0] * (shape[1] - 1)


def vertical_link_count(shape):
    """Number of vertical links."""
    assert len(shape) == 2
    return (shape[0] - 1) * shape[1]


def horizontal_link_count(shape):
    """Number of horizontal links."""
    assert len(shape) == 2
    return shape[0] * (shape[1] - 1)


def perimeter_node_count(shape):
    """Number of perimeter nodes.

    Number of nodes that are on the perimeter of a structured grid with
    dimensions, *shape*, and thus boundary nodes.

    Examples
    --------
    >>> from landlab.utils.structured_grid import perimeter_node_count
    >>> perimeter_node_count((3, 4))
    10
    """
    assert len(shape) == 2
    return 2 * (shape[0] - 2) + 2 * (shape[1] - 2) + 4


def interior_cell_count(shape):
    """Number of interior cells.

    Number of interior cells. Since cells are only defined on interior nodes,
    this is the same as cell_count.
    """
    return cell_count(shape)


def face_count(shape):
    """Total number of faces.

    Total number of faces in a structured grid with dimensions, *shape*. Each
    cell has four faces, and shared faces only count once.

    Examples
    --------
    >>> from landlab.utils.structured_grid import face_count
    >>> face_count((3, 4))
    7
    """
    assert len(shape) == 2
    if np.min(shape) > 2:
        return (shape[0] - 1) * (shape[1] - 2) + (shape[0] - 2) * (shape[1] - 1)
    else:
        return 0


def active_face_count(shape):
    """Number of active faces.

    Total number of active faces in a structured grid with dimensions,
    *shape*. Each cell has four faces, and shared faces only count once.
    An active face is one that has a corresponing active link.

    >>> from landlab.utils.structured_grid import active_face_count
    >>> active_face_count((3, 4))
    7
    """
    return face_count(shape)


def top_index_iter(shape):
    """Iterator for the top boundary indices of a structured grid."""
    return range(shape[1] * (shape[0] - 1), shape[0] * shape[1])


def bottom_index_iter(shape):
    """Iterator for the bottom boundary indices of a structured grid."""
    return range(0, shape[1])


def left_index_iter(shape):
    """Iterator for the left boundary indices of a structured grid."""
    return range(0, shape[0] * shape[1], shape[1])


def right_index_iter(shape):
    """Iterator for the right boundary indices of a structured grid."""
    return range(shape[1] - 1, shape[0] * shape[1], shape[1])


def left_right_iter(shape, *args):
    """Iterator for the left and right boundary indices of a structured grid.

    This iterates over the indices in order rather than iterating all of
    the left boundary and then all of the right boundary.

    Examples
    --------
    >>> import numpy as np
    >>> from landlab.utils.structured_grid import left_right_iter
    >>> np.fromiter(left_right_iter((4, 3)), dtype=np.int)
    array([ 0,  2,  3,  5,  6,  8,  9, 11])
    >>> np.fromiter(left_right_iter((4, 3), 2), dtype=np.int)
    array([0, 2, 3, 5])
    >>> np.fromiter(left_right_iter((4, 3), 2, 4), dtype=np.int)
    array([ 6,  8,  9, 11])
    >>> np.fromiter(left_right_iter((4, 3), 1, 4, 2), dtype=np.int)
    array([ 3,  5,  9, 11])
    """
    if len(args) == 0:
        iter_rows = range(0, shape[0], 1)
    elif len(args) == 1:
        iter_rows = range(0, args[0], 1)
    elif len(args) == 2:
        iter_rows = range(args[0], args[1], 1)
    elif len(args) == 3:
        iter_rows = range(args[0], args[1], args[2])

    for row in iter_rows:
        yield row * shape[1]
        yield row * shape[1] + shape[1] - 1


def bottom_top_iter(shape):
    """Iterator for bottom then top indices of a structured grid.

    Examples
    --------
    >>> import numpy as np
    >>> from landlab.utils.structured_grid import bottom_top_iter
    >>> np.fromiter(bottom_top_iter((4, 3)), dtype=np.int)
    array([ 0,  1,  2,  9, 10, 11])
    """
    return itertools.chain(bottom_index_iter(shape), top_index_iter(shape))


def boundary_iter(shape):
    """Iterator for perimeter nodes.

    .. deprecated:: 0.6

        Deprecated due to imprecise terminology. This is really perimeter_iter
        (see below).
        Iterates over all of the boundary node indices of a structured grid in
        order.

    Examples
    --------
    >>> import numpy as np
    >>> from landlab.utils.structured_grid import boundary_iter
    >>> np.fromiter(boundary_iter((4, 3)), dtype=np.int)
    array([ 0,  1,  2,  3,  5,  6,  8,  9, 10, 11])
    """
    return itertools.chain(
        bottom_index_iter(shape),
        left_right_iter(shape, 1, shape[0] - 1),
        top_index_iter(shape),
    )


def perimeter_iter(shape):
    """Iterator for perimeter nodes.

    Iterates over all of the perimeter node indices of a structured grid in
    order.

    Examples
    --------
    >>> import numpy as np
    >>> from landlab.utils.structured_grid import perimeter_iter
    >>> np.fromiter(perimeter_iter((4, 3)), dtype=np.int)
    array([ 0,  1,  2,  3,  5,  6,  8,  9, 10, 11])
    """
    return itertools.chain(
        bottom_index_iter(shape),
        left_right_iter(shape, 1, shape[0] - 1),
        top_index_iter(shape),
    )


def boundary_nodes(shape):
    """Array of perimeter nodes.

    .. deprecated:: 0.6
        Deprecated due to imprecise terminology. This is really perimeter_iter
        (see below).

    An array of the indices of the boundary nodes.

    Examples
    --------
    >>> from landlab.utils.structured_grid import boundary_nodes
    >>> boundary_nodes((3, 4))
    array([ 0,  1,  2,  3,  4,  7,  8,  9, 10, 11])
    """
    return np.fromiter(boundary_iter(shape), dtype=np.int)


def perimeter_nodes(shape):
    """Array of perimeter nodes.

    An array of the indices of the perimeter nodes of a structured grid.

    Examples
    --------
    >>> from landlab.utils.structured_grid import perimeter_nodes
    >>> perimeter_nodes((3, 4))
    array([ 0,  1,  2,  3,  4,  7,  8,  9, 10, 11])
    """
    return np.fromiter(perimeter_iter(shape), dtype=np.int)


def corners(shape):
    """Array of the indices of the grid corner nodes."""
    return np.array(
        [0, shape[1] - 1, shape[1] * (shape[0] - 1), shape[1] * shape[0] - 1]
    )


def bottom_edge_node_ids(shape):
    """Array of nodes on the bottom edge."""
    return np.fromiter(bottom_index_iter(shape), dtype=np.int)


def top_edge_node_ids(shape):
    """Array of nodes on the top edge."""
    return np.fromiter(top_index_iter(shape), dtype=np.int)


def left_edge_node_ids(shape):
    """Array of nodes on the left edge."""
    return np.fromiter(left_index_iter(shape), dtype=np.int)


def right_edge_node_ids(shape):
    """Array of nodes on the right edge."""
    return np.fromiter(right_index_iter(shape), dtype=np.int)


def interior_iter(shape):
    """Iterator for the interior nodes of a structured grid.

    Examples
    --------
    >>> import numpy as np
    >>> from landlab.utils.structured_grid import interior_iter
    >>> np.fromiter(interior_iter((4, 3)), dtype=np.int)
    array([4, 7])
    """
    interiors = []
    interiors_per_row = shape[1] - 2
    for row in range(shape[1] + 1, shape[1] * (shape[0] - 1), shape[1]):
        interiors.append(range(row, row + interiors_per_row))
    return itertools.chain(*interiors)


def interior_nodes(shape):
    """Array of interior nodes."""
    return np.fromiter(interior_iter(shape), dtype=np.int)


def node_coords(shape, *args):
    """Get node x and y coordinates.

    Get x, y coordinates for nodes in a structured grid with dimensions,
    *shape*. Use the optional argument *spacing* to give the spacing in each
    dimension, and *origin* the start of the coordinates in each dimension.

    Parameters
    ----------
    shape : tuple of int
        Number of node rows and columns.
    spacing : tuple, optional
        Row and column spacing.
    origin : tuple, optional
        Coordinate of lower-left node.

    Examples
    --------
    >>> from landlab.utils.structured_grid import node_coords
    >>> (cols, rows) = node_coords((3, 2))
    >>> rows
    array([ 0.,  0.,  1.,  1.,  2.,  2.])
    >>> cols
    array([ 0.,  1.,  0.,  1.,  0.,  1.])
    """
    try:
        spacing = args[0]
    except IndexError:
        spacing = np.ones(len(shape), dtype=np.float)
    else:
        assert len(spacing) == len(shape)

    try:
        origin = args[1]
    except IndexError:
        origin = np.zeros(len(shape), dtype=np.float)
    else:
        assert len(origin) == len(origin)

    node_count_ = np.prod(shape)

    row_y = np.arange(shape[0]) * spacing[0] + origin[0]
    col_x = np.arange(shape[1]) * spacing[1] + origin[1]

    (node_x, node_y) = np.meshgrid(col_x, row_y)

    node_x.shape = (node_count_,)
    node_y.shape = (node_count_,)

    return (node_x, node_y)


def node_at_cell(shape):
    """Array of nodes at cells.

    Indices of the nodes belonging to each cell.

    Examples
    --------
    >>> from landlab.utils.structured_grid import node_at_cell
    >>> node_at_cell((4, 3))
    array([4, 7])
    """
    node_ids = np.arange(node_count(shape))
    node_ids.shape = shape

    cell_node = node_ids[1:-1, 1:-1].copy()
    cell_node.shape = ((shape[0] - 2) * (shape[1] - 2),)

    return cell_node


def node_index_at_link_ends(shape):
    """Array of nodes at each end of links."""
    node_ids = np.arange(np.prod(shape))
    node_ids.shape = shape

    return (node_at_link_tail(node_ids), node_at_link_head(node_ids))


def inlink_index_at_node(shape):
    """Array of links entering nodes."""
    return inlinks(shape)


def outlink_index_at_node(shape):
    """Array of links leaving nodes."""
    return outlinks(shape)


def node_at_link_head(node_ids):
    """Array of nodes at the end of links."""
    vertical_links = node_ids[1:, :]
    horizontal_links = node_ids[:, 1:]
    return np.concatenate((vertical_links.flat, horizontal_links.flat))


def node_at_link_tail(node_ids):
    """Array of nodes at the start of links."""
    vertical_links = node_ids[:-1, :]
    horizontal_links = node_ids[:, :-1]
    return np.concatenate((vertical_links.flat, horizontal_links.flat))


def face_at_link(shape, actives=None, inactive_link_index=BAD_INDEX_VALUE):
    """Array of faces associated with links.

    Returns an array that maps link ids to face ids. For inactive links,
    which do not have associated faces, set their ids to
    *inactive_link_index*. Use the *actives* keyword to specify an array that
    contains the ids of all active links in the grid. The default assumes
    that only the perimeter nodes are inactive.

    Examples
    --------
    >>> from landlab.utils.structured_grid import face_at_link
    >>> faces = face_at_link((3, 4), inactive_link_index=-1)
    >>> faces # doctest: +NORMALIZE_WHITESPACE
    array([-1,  0,  1, -1, -1,  2,  3,
           -1, -1, -1, -1,  4,  5,  6, -1, -1, -1])
    """
    if actives is None:
        actives = active_links(shape)

    num_links = link_count(shape)

    link_faces = np.empty(num_links, dtype=np.int)
    link_faces.fill(inactive_link_index)
    link_faces[actives] = np.arange(len(actives))

    return link_faces


def status_at_node(shape, boundary_status=FIXED_VALUE_BOUNDARY):
    """Array of the statuses of nodes.

    The statuses of the nodes in a structured grid with dimensions, *shape*.
    Use the *boundary_status* keyword to specify the status of the top,
    bottom, left and right boundary nodes.
    """
    status = np.empty(np.prod(shape), dtype=np.int8)

    status[interior_nodes(shape)] = CORE_NODE
    status[boundary_nodes(shape)] = boundary_status

    return status


def active_links(shape, node_status_array=None, link_nodes=None):
    """Link IDs for active links of a structured quad grid.

    Return the link IDs for links that are *active* in a structured grid of
    quadrilaterals. Use the *node_status_array* keyword to specify the status
    for each of the grid's nodes. If not given, each of the perimeter nodes is
    assumed to be `FIXED_VALUE_BOUNDARY`.

    Use the *link_nodes* keyword to provide, as a tuple of arrays, that give
    the *from-node* and the *to-node* for each for each link in the grid.

    Parameters
    ----------
    shape : tuple
        Shape of grid as number of node rows and columns.
    node_status_array : array_like, optional
        Status of each grid node.
    link_nodes : array_like, optional

    Examples
    --------
    Because, by default, the perimeter nodes are `FIXED_VALUE_BOUNDARY` nodes,
    only links attached to the interior nodes are *active*.

    >>> from landlab.utils.structured_grid import active_links
    >>> from landlab import CLOSED_BOUNDARY, CORE_NODE
    >>> active_links((3, 4))
    array([ 1,  2,  5,  6, 11, 12, 13])

    If all the perimeter nodes `CLOSED_BOUNDARY` nodes, the only active link
    is between the two core nodes.

    >>> node_status = np.ones(3 * 4) * CLOSED_BOUNDARY
    >>> node_status[5:7] = CORE_NODE
    >>> active_links((3, 4), node_status_array=node_status)
    array([12])

    You can also provide a list of all the *from_nodes* and *to_nodes* for
    the grid. The following describes a grid with only a single link (between
    nodes 5 and 6).

    >>> active_links((3, 4), link_nodes=(np.array([5]), np.array([6])))
    array([0])
    """
    if node_status_array is None:
        node_status_array = status_at_node(shape)

    if link_nodes is None:
        (link_from_node, link_to_node) = node_index_at_link_ends(shape)
    else:
        (link_from_node, link_to_node) = link_nodes

    from_node_status = node_status_array[link_from_node]
    to_node_status = node_status_array[link_to_node]

    active_links_ = (
        (from_node_status == CORE_NODE) & ~(to_node_status == CLOSED_BOUNDARY)
    ) | ((to_node_status == CORE_NODE) & ~(from_node_status == CLOSED_BOUNDARY))

    (active_links_,) = np.where(active_links_)

    return as_id_array(active_links_)


def active_face_index(shape):
    """Array of face indices."""
    return np.arange(active_face_count(shape))


def inlinks(shape):
    """Array of links entering nodes."""
    links = np.vstack((south_links(shape), west_links(shape)))
    links.shape = (2, node_count(shape))
    return links


def outlinks(shape):
    """Array of links leaving nodes."""
    links = np.vstack((north_links(shape), east_links(shape)))
    links.shape = (2, node_count(shape))
    return links


def active_inlinks(shape, node_status=None):
    """Array of active links entering nodes."""
    links = np.vstack(
        (active_south_links(shape, node_status), active_west_links(shape, node_status))
    )
    links.shape = (2, node_count(shape))
    return links


def active_inlinks2(shape, node_status=None):
    """Array of active links entering nodes.

    Finds and returns the link IDs of active links coming in to each node
    (that is, active links for which the node is the link head).

    Parameters
    ----------
    shape : 2-element tuple of ints
        Number of rows and columns in the grid
    node_status (optional) : numpy array of bool (x # of nodes)
        False where node is a closed boundary; True elsewhere

    Returns
    -------
    2d numpy array of int (2 x number of grid nodes)
        Link ID of incoming links to each node

    Examples
    --------
    >>> from landlab.utils.structured_grid import active_inlinks2
    >>> active_inlinks2((3,4))
    array([[-1, -1, -1, -1, -1,  4,  5, -1, -1, 11, 12, -1],
           [-1, -1, -1, -1, -1,  7,  8,  9, -1, -1, -1, -1]])

    Notes
    -----
    There are at most two inlinks for each node. The first row in the returned
    array gives the ID of the vertical incoming link from below (south), or -1
    if there is none. The second row gives the link ID of the horizontal link
    coming in from the left (or -1).
    """
    links = np.vstack(
        (
            active_south_links2(shape, node_status),
            active_west_links2(shape, node_status),
        )
    )
    links.shape = (2, node_count(shape))
    return links


def active_outlinks(shape, node_status=None):
    """Array of active links leaving nodes."""
    links = np.vstack(
        (active_north_links(shape, node_status), active_east_links(shape, node_status))
    )
    links.shape = (2, node_count(shape))
    return links


def active_outlinks2(shape, node_status=None):
    """Array of active links leaving nodes.

    Finds and returns the link IDs of active links going out of each node
    (that is, active links for which the node is the link tail).

    Parameters
    ----------
    shape : 2-element tuple of ints
        Number of rows and columns in the grid
    node_status (optional) : numpy array of bool (x # of nodes)
        False where node is a closed boundary; True elsewhere

    Returns
    -------
    2d numpy array of int (2 x number of grid nodes)
        Link ID of outgoing links from each node

    Examples
    --------
    >>> from landlab.utils.structured_grid import active_outlinks2
    >>> active_outlinks2((3,4))
    array([[-1,  4,  5, -1, -1, 11, 12, -1, -1, -1, -1, -1],
           [-1, -1, -1, -1,  7,  8,  9, -1, -1, -1, -1, -1]])

    Notes
    -----
    There are at most two outlinks for each node. The first row in the returned
    array gives the ID of the vertical outgoing link to above (north), or -1
    if there is none. The second row gives the link ID of the horizontal link
    going out to the right (east) (or -1).
    """
    links = np.vstack(
        (
            active_north_links2(shape, node_status),
            active_east_links2(shape, node_status),
        )
    )
    links.shape = (2, node_count(shape))
    return links


def vertical_link_ids(shape):
    """Array of links oriented vertically."""
    link_ids = np.empty((shape[0] - 1, shape[1]), dtype=np.int)
    num_links_per_row = (2 * shape[1]) - 1
    for r in range(shape[0] - 1):
        link_ids[r, :] = (shape[1] - 1) + (r * num_links_per_row) + np.arange(shape[1])
    return link_ids


def horizontal_link_ids(shape):
    """Array of links oriented horizontally."""
    link_ids = np.empty((shape[0], shape[1] - 1), dtype=np.int)
    num_links_per_row = (2 * shape[1]) - 1
    for r in range(shape[0]):
        link_ids[r, :] = (r * num_links_per_row) + np.arange(shape[1] - 1)
    return link_ids


def vertical_active_link_count(shape, node_status=None):
    """Number of active links oriented vertically."""
    if node_status is not None:
        is_inactive = node_status == 0
        is_inactive.shape = shape

        inactive_outlinks = is_inactive[:-1, 1:-1]
        inactive_inlinks = is_inactive[1:, 1:-1]
        inactive_links = inactive_outlinks | inactive_inlinks

        return (shape[0] - 1) * (shape[1] - 2) - np.sum(inactive_links.flat)
    else:
        return (shape[0] - 1) * (shape[1] - 2)


def horizontal_active_link_count(shape, node_status=None):
    """Number of active links oriented horizontally."""
    if node_status is not None:
        is_inactive = node_status == 0
        is_inactive.shape = shape

        inactive_outlinks = is_inactive[1:-1, :-1]
        inactive_inlinks = is_inactive[1:-1, 1:]
        inactive_links = inactive_outlinks | inactive_inlinks

        return (shape[0] - 2) * (shape[1] - 1) - np.sum(inactive_links.flat)
    else:
        return (shape[0] - 2) * (shape[1] - 1)


def vertical_inactive_link_mask(shape, node_status):
    """Array mask of vertical links that are inactive.

    Creates and returns a boolean 2D array dimensioned as the number of
    vertical links in the grid, not including the left and right boundaries.

    Parameters
    ----------
    shape : 2-element tuple of ints
        Number of rows and columns in the grid
    node_status : numpy array of bool (x # of nodes)
        False where node is a closed boundary; True elsewhere

    Returns
    -------
    (NR-1,NC-2) array of bool (NR=# of rows, NC=# of columns)
        Flags indicating whether the corresponding vertical link is inactive

    Examples
    --------
    >>> import numpy as np
    >>> from landlab.utils.structured_grid import vertical_inactive_link_mask
    >>> ns = np.ones(12, dtype=bool)  # case of no closed boundary nodes
    >>> vertical_inactive_link_mask((3,4), ns)
    array([[False, False],
           [False, False]], dtype=bool)
    >>> ns[2] = False    # node 2 is a closed boundary
    >>> vertical_inactive_link_mask((3,4), ns)
    array([[False,  True],
           [False, False]], dtype=bool)
    >>> ns[9] = False    # node 9 is also a closed boundary
    >>> vertical_inactive_link_mask((3,4), ns)
    array([[False,  True],
           [ True, False]], dtype=bool)
    """
    # Create a 2D boolean matrix indicating whether NODES are closed boundaries
    # GT thinks this should be False, not 0
    is_closed_node = node_status == 0
    is_closed_node.shape = shape

    inactive_outlinks = is_closed_node[:-1, 1:-1]  # middle cols, all but top row
    # middle cols, all but bottom row
    inactive_inlinks = is_closed_node[1:, 1:-1]
    # if either node is closed, the link is inactive
    return inactive_outlinks | inactive_inlinks


def horizontal_inactive_link_mask(shape, node_status):
    """Array mask of horizontal links that are inactive."""
    is_inactive = node_status == 0
    is_inactive.shape = shape

    inactive_outlinks = is_inactive[1:-1, :-1]
    inactive_inlinks = is_inactive[1:-1, 1:]
    return inactive_outlinks | inactive_inlinks


# def vertical_active_link_ids(shape):
#    link_ids = np.arange(vertical_active_link_count(shape), dtype=np.int)
#    link_ids.shape = (shape[0] - 1, shape[1] - 2)
#    return link_ids


def vertical_active_link_ids(shape, node_status=None):
    """Array of active links oriented vertically."""
    if node_status is None:
        link_ids = np.arange(vertical_active_link_count(shape), dtype=np.int)
        # link_ids.shape = (shape[0] - 1, shape[1] - 2)
        # return link_ids
    else:
        inactive_links = vertical_inactive_link_mask(shape, node_status)
        inactive_links.shape = (inactive_links.size,)
        active_link_count_ = inactive_links.size - np.sum(inactive_links)

        link_ids = np.empty(inactive_links.size)
        link_ids[inactive_links] = -1
        link_ids[~inactive_links] = np.arange(active_link_count_, dtype=int)

    link_ids.shape = (shape[0] - 1, shape[1] - 2)
    return link_ids


def vertical_active_link_ids2(shape, node_status=None):
    """Array of active links oriented vertically.

    Returns the link IDs of vertical active links as an (R-1) x (C-2) array.

    Parameters
    ----------
    shape : 2-element tuple of int
        number of rows and columns in grid
    node_status (optional) : 1d numpy array (x number of nodes) of bool
        False where node is a closed boundary, True otherwise

    Returns
    -------
    2d numpy array of int
        Link IDs of vertical active links, not including vertical links on the
        left and right grid edges. If a vertical link is inactive, its ID is
        given as -1.

    Examples
    --------
    >>> from landlab.utils.structured_grid import vertical_active_link_ids2
    >>> vertical_active_link_ids2((3,4))
    array([[ 4,  5],
           [11, 12]])
    >>> ns = np.ones(12, dtype=bool)
    >>> ns[1] = False
    >>> ns[10] = False
    >>> vertical_active_link_ids2((3,4), ns)
    array([[-1,  5],
           [11, -1]])

    Notes
    -----
    Same as vertical_active_link_ids() but returns "link IDs" for active links
    rather than "active link IDs" for active links. Designed to ultimately
    replace the original vertical_active_link_ids().
    """
    link_ids = np.empty((shape[0] - 1, shape[1] - 2), dtype=np.int)
    num_links_per_row = (2 * shape[1]) - 1
    for r in range(shape[0] - 1):
        link_ids[r, :] = shape[1] + (r * num_links_per_row) + np.arange(shape[1] - 2)

    if node_status is not None:
        inactive_links = vertical_inactive_link_mask(shape, node_status)
        link_ids[inactive_links] = -1

    return link_ids


def horizontal_active_link_ids(shape, node_status=None):
    """Array of active links oriented horizontally."""
    if node_status is None:
        link_id_offset = vertical_active_link_count(shape)
        link_ids = np.arange(
            link_id_offset,
            link_id_offset + horizontal_active_link_count(shape),
            dtype=np.int,
        )
    else:
        link_id_offset = vertical_active_link_count(shape, node_status=node_status)
        inactive_links = horizontal_inactive_link_mask(shape, node_status)
        inactive_links.shape = (inactive_links.size,)
        active_link_count_ = inactive_links.size - np.sum(inactive_links)

        link_ids = np.empty(inactive_links.size)
        link_ids[inactive_links] = -1
        link_ids[~inactive_links] = np.arange(
            link_id_offset, link_id_offset + active_link_count_, dtype=np.int
        )
    link_ids.shape = (shape[0] - 2, shape[1] - 1)
    return link_ids


def horizontal_active_link_ids2(shape, node_status=None):
    """Array of active links oriented horizontally.

    Returns the link IDs of horizontal active links as an (R-2) x (C-1) array.

    Parameters
    ----------
    shape : 2-element tuple of int
        number of rows and columns in grid
    node_status (optional) : 1d numpy array (x number of nodes) of bool
        False where node is a closed boundary, True otherwise

    Returns
    -------
    2d numpy array of int
        Link IDs of horizontal active links, not including horizontal links on
        top and bottom grid edges. If a horizontal link is inactive, its ID is
        given as -1.

    Examples
    --------
    >>> from landlab.utils.structured_grid import horizontal_active_link_ids2
    >>> horizontal_active_link_ids2((3,4))
    array([[7, 8, 9]])
    >>> ns = np.ones(12, dtype=bool)
    >>> ns[4] = False
    >>> ns[7] = False
    >>> horizontal_active_link_ids2((3,4), ns)
    array([[-1,  8, -1]])

    Notes
    -----
    Same as horizontal_active_link_ids() but returns "link IDs" for active
    links rather than "active link IDs" for active links. Designed to
    ultimately replace the original horizontal_active_link_ids().
    """
    link_ids = np.empty((shape[0] - 2, shape[1] - 1), dtype=np.int)
    num_links_per_row = (2 * shape[1]) - 1
    for r in range(shape[0] - 2):
        link_ids[r, :] = ((r + 1) * num_links_per_row) + np.arange(shape[1] - 1)
    if node_status is not None:
        inactive_links = horizontal_inactive_link_mask(shape, node_status)
        link_ids[inactive_links] = -1

    return link_ids


def west_links(shape):
    """Array of links pointing to the west.

    Examples
    --------
    >>> from landlab.utils.structured_grid import west_links
    >>> west_links((3, 4))
    array([[-1,  0,  1,  2],
           [-1,  7,  8,  9],
           [-1, 14, 15, 16]])
    """
    link_ids = horizontal_link_ids(shape)
    link_ids.shape = (shape[0], shape[1] - 1)
    return np.hstack((-np.ones((shape[0], 1), dtype=np.int), link_ids))


def north_links(shape):
    """Array of links pointing to the north.

    Examples
    --------
    >>> from landlab.utils.structured_grid import north_links
    >>> north_links((3, 4))
    array([[ 3,  4,  5,  6],
           [10, 11, 12, 13],
           [-1, -1, -1, -1]])
    """
    link_ids = vertical_link_ids(shape)
    link_ids.shape = (shape[0] - 1, shape[1])
    return np.vstack((link_ids, -np.ones((1, shape[1]), dtype=np.int)))


def south_links(shape):
    """Array of links pointing to the the south.

    Examples
    --------
    >>> from landlab.utils.structured_grid import south_links
    >>> south_links((3, 4))
    array([[-1, -1, -1, -1],
           [ 3,  4,  5,  6],
           [10, 11, 12, 13]])
    """
    link_ids = vertical_link_ids(shape)
    link_ids.shape = (shape[0] - 1, shape[1])
    return np.vstack((-np.ones((1, shape[1]), dtype=np.int), link_ids))


def east_links(shape):
    """Array of links pointing to the the east.

    Examples
    --------
    >>> from landlab.utils.structured_grid import east_links
    >>> east_links((3, 4))
    array([[ 0,  1,  2, -1],
           [ 7,  8,  9, -1],
           [14, 15, 16, -1]])
    """
    link_ids = horizontal_link_ids(shape)
    link_ids.shape = (shape[0], shape[1] - 1)
    return np.hstack((link_ids, -np.ones((shape[0], 1), dtype=int)))


def active_north_links(shape, node_status=None):
    """Array of active links pointing to the the north."""
    active_north_links_ = np.empty(shape, dtype=int)
    try:
        links = vertical_active_link_ids(shape, node_status=node_status)
    except ValueError:
        pass
    links.shape = (shape[0] - 1, shape[1] - 2)
    active_north_links_[:-1, 1:-1] = links
    active_north_links_[:, (0, -1)] = -1
    active_north_links_[-1, :] = -1

    return active_north_links_


def active_north_links2(shape, node_status=None):
    """Array of active links pointing to the the north.

    >>> from landlab.utils.structured_grid import active_north_links2
    >>> active_north_links2((3, 4))
    array([[-1,  4,  5, -1],
           [-1, 11, 12, -1],
           [-1, -1, -1, -1]])
    """
    active_north_links_ = np.empty(shape, dtype=int)
    try:
        links = vertical_active_link_ids2(shape, node_status=node_status)
    except ValueError:
        pass
    links.shape = (shape[0] - 1, shape[1] - 2)
    active_north_links_[:-1, 1:-1] = links
    active_north_links_[:, (0, -1)] = -1
    active_north_links_[-1, :] = -1

    return active_north_links_


def active_south_links(shape, node_status=None):
    """Array of active links pointing to the the south."""
    active_south_links_ = np.empty(shape, dtype=int)
    links = vertical_active_link_ids(shape, node_status=node_status)
    links.shape = (shape[0] - 1, shape[1] - 2)
    active_south_links_[1:, 1:-1] = links
    active_south_links_[:, (0, -1)] = -1
    active_south_links_[0, :] = -1

    return active_south_links_


def active_south_links2(shape, node_status=None):
    """Array of active links pointing to the the south.

    Finds and returns link IDs of active links that enter each node from the
    south (bottom), or -1 where no such active link exists.

    Parameters
    ----------
    shape : 2-element tuple of int
        number of rows and columns in grid
    node_status (optional) : 1d numpy array of bool
        False where node is a closed boundary, True otherwise

    Returns
    -------
    2d numpy array of int
        Link ID of active link connecting to a node from the south, or -1

    Examples
    --------
    >>> from landlab.utils.structured_grid import active_south_links2
    >>> active_south_links2((3, 4))
    array([[-1, -1, -1, -1],
           [-1,  4,  5, -1],
           [-1, 11, 12, -1]])

    Notes
    -----
    Like active_south_links, but returns link IDs.
    """
    active_south_links_ = -np.ones(shape, dtype=int)
    links = vertical_active_link_ids2(shape, node_status=node_status)
    active_south_links_[1:, 1:-1] = links

    return active_south_links_


def active_west_links(shape, node_status=None):
    """Array of active links pointing to the the west.

    Examples
    --------
    >>> from landlab.utils.structured_grid import active_west_links
    >>> active_west_links((3, 4))
    array([[-1, -1, -1, -1],
           [-1,  4,  5,  6],
           [-1, -1, -1, -1]])
    """
    active_west_links_ = np.empty(shape, dtype=int)
    try:
        active_west_links_[1:-1, 1:] = horizontal_active_link_ids(
            shape, node_status=node_status
        )
    except ValueError:
        pass
    active_west_links_[(0, -1), :] = -1
    active_west_links_[:, 0] = -1

    return active_west_links_


def active_west_links2(shape, node_status=None):
    """Array of active links pointing to the the west.

    Examples
    --------
    >>> from landlab.utils.structured_grid import active_west_links2
    >>> active_west_links2((3, 4))
    array([[-1, -1, -1, -1],
           [-1,  7,  8,  9],
           [-1, -1, -1, -1]])
    """
    active_west_links_ = -np.ones(shape, dtype=int)
    active_west_links_[1:-1, 1:] = horizontal_active_link_ids2(
        shape, node_status=node_status
    )

    return active_west_links_


def active_east_links(shape, node_status=None):
    """Array of active links pointing to the the east.

    Examples
    --------
    >>> from landlab.utils.structured_grid import active_east_links
    >>> active_east_links((3, 4))
    array([[-1, -1, -1, -1],
           [ 4,  5,  6, -1],
           [-1, -1, -1, -1]])
    """
    active_east_links_ = np.empty(shape, dtype=int)
    active_east_links_.fill(-999)
    try:
        active_east_links_[1:-1, :-1] = horizontal_active_link_ids(
            shape, node_status=node_status
        )
    except ValueError:
        pass
    active_east_links_[(0, -1), :] = -1
    active_east_links_[:, -1] = -1

    return active_east_links_


def active_east_links2(shape, node_status=None):
    """Array of active links pointing to the the east.

    Examples
    --------
    >>> from landlab.utils.structured_grid import active_east_links2
    >>> active_east_links2((3, 4))
    array([[-1, -1, -1, -1],
           [ 7,  8,  9, -1],
           [-1, -1, -1, -1]])
    """
    active_east_links_ = -np.ones(shape, dtype=int)
    active_east_links_[1:-1, :-1] = horizontal_active_link_ids2(
        shape, node_status=node_status
    )

    return active_east_links_


def outlink_count_per_node(shape):
    """Number of links leaving each node."""
    link_count_ = np.empty(shape, dtype=np.int)
    link_count_[:-1, :-1] = 2
    link_count_[-1, :-1] = 1
    link_count_[:-1, -1] = 1
    link_count_[-1, -1] = 0
    return np.ravel(link_count_)


def inlink_count_per_node(shape):
    """Number of links entering each node."""
    link_count_ = np.empty(shape, dtype=np.int)
    link_count_[1:, 1:] = 2
    link_count_[0, 1:] = 1
    link_count_[1:, 0] = 1
    link_count_[0, 0] = 0
    return np.ravel(link_count_)


def active_outlink_count_per_node(shape):
    """Number of active links leaving each node."""
    link_count_ = np.empty(shape, dtype=np.int)
    link_count_[1:-1, 1:-1] = 2
    link_count_[0, :] = 1
    link_count_[-1, :] = 0
    link_count_[:, 0] = 1
    link_count_[:, -1] = 0

    link_count_[0, 0] = 0
    link_count_[-1, 0] = 0

    return np.ravel(link_count_)


def active_inlink_count_per_node(shape):
    """Number of active links entering each node."""
    link_count_ = np.empty(shape, dtype=np.int)
    link_count_[1:-1, 1:-1] = 2
    link_count_[0, :] = 0
    link_count_[-1, :] = 1
    link_count_[:, 0] = 0
    link_count_[:, -1] = 1

    link_count_[0, -1] = 0
    link_count_[-1, -1] = 0

    return np.ravel(link_count_)


def setup_outlink_matrix(shape, return_count=True):
    """Create a matrix of links leaving each node."""
    links = outlinks(shape)
    if return_count:
        return (links, outlink_count_per_node(shape))
    else:
        return links


def setup_inlink_matrix(shape, return_count=True):
    """Create a matrix of links entering each node."""
    links = inlinks(shape)
    if return_count:
        return (links, inlink_count_per_node(shape))
    else:
        return links


def setup_active_outlink_matrix(shape, node_status=None, return_count=True):
    """Create a matrix of active links leaving each node."""
    links = active_outlinks(shape, node_status=node_status)
    if return_count:
        return links, active_outlink_count_per_node(shape)
    else:
        return links


def setup_active_outlink_matrix2(shape, node_status=None, return_count=True):
    """Create a matrix of active links leaving each node.

    Return the link IDs of the active links that leave each node of a grid. The
    shape of the returned array is (2, *N*) where *N* is the number of nodes in
    the grid. The first row contains the link ID exiting the node to the
    top, and the second row the link exiting the node to the right.

    Use the *return_count* keyword to, in addition to the link IDs, return the
    number of active links attached to each grid node.

    Use the *node_status_array* keyword to specify the status for each of the
    grid's nodes. If not given, each of the perimeter nodes is assumed to be
    `FIXED_VALUE_BOUNDARY`.

    Parameters
    ----------
    shape : tuple
        Shape of the structured grid
    node_status : array_like, optional
        Status of each node in the grid.
    return_count : boolean, optional
        If `True`, also return an array of active link counts per node.

    Returns
    -------
    links : (2, N) ndarray
        Active link IDs for each node.
    count : ndarray
        Number of active links per node.

    Examples
    --------
    Get the active link IDs for a grid of 3 nodes by 4 nodes. The first row
    lists links entering nodes from the bottom, and the second links entering
    from the left.

    >>> from landlab.utils.structured_grid import setup_active_outlink_matrix2
    >>> setup_active_outlink_matrix2((3, 4), return_count=False)
    array([[-1,  4,  5, -1, -1, 11, 12, -1, -1, -1, -1, -1],
           [-1, -1, -1, -1,  7,  8,  9, -1, -1, -1, -1, -1]])
    >>> _, count = setup_active_outlink_matrix2((3, 4))
    >>> count
    array([0, 1, 1, 0, 1, 2, 2, 0, 0, 0, 0, 0])
    """
    links = active_outlinks2(shape, node_status=node_status)
    if return_count:
        return links, active_outlink_count_per_node(shape)
    else:
        return links


def setup_active_inlink_matrix(shape, node_status=None, return_count=True):
    """Create a matrix of active links entering each node.

    Return the IDs of the active links that enter each node of a grid. The
    shape of the returned array is (2, *N*) where *N* is the number of nodes
    in the grid. The first row contains the link ID entering the node from the
    bottom, and the second row the link entering the node from the left.

    Use the *return_count* keyword to, in addition to the link IDs, return
    the number of active links attached to each grid node.

    Use the *node_status_array* keyword to specify the status for each of
    the grid's nodes. If not given, each of the perimeter nodes is assumed
    to be `FIXED_VALUE_BOUNDARY`.

    Parameters
    ----------
    shape : tuple
        Shape of the structured grid
    node_status : array_like, optional
        Status of each node in the grid.
    return_count : boolean, optional
        If `True`, also return an array of active link counts per node.

    Returns
    -------
    links : (2, N) ndarray
        Active link IDs for each node.
    count : ndarray
        Number of active links per node.

    Examples
    --------
    Get the active link IDs for a grid of 3 nodes by 4 nodes. The first row
    list links entering nodes from the bottom, and the second links entering
    from the left.

    >>> from landlab.utils.structured_grid import setup_active_inlink_matrix
    >>> setup_active_inlink_matrix((3, 4), return_count=False)
    array([[-1, -1, -1, -1, -1,  0,  1, -1, -1,  2,  3, -1],
           [-1, -1, -1, -1, -1,  4,  5,  6, -1, -1, -1, -1]])
    >>> _, count = setup_active_inlink_matrix((3, 4))
    >>> count
    array([0, 0, 0, 0, 0, 2, 2, 1, 0, 1, 1, 0])
    """
    links = active_inlinks(shape, node_status=node_status)
    if return_count:
        return links, active_inlink_count_per_node(shape)
    else:
        return links


def setup_active_inlink_matrix2(shape, node_status=None, return_count=True):
    """Create a matrix of active links entering each node.

    Return the link IDs of the active links that enter each node of a grid. The
    shape of the returned array is (2, *N*) where *N* is the number of nodes in
    the grid. The first row contains the link ID entering the node from the
    bottom, and the second row the link entering the node from the left.

    Use the *return_count* keyword to, in addition to the link IDs, return the
    number of active links attached to each grid node.

    Use the *node_status_array* keyword to specify the status for each of the
    grid's nodes. If not given, each of the perimeter nodes is assumed to be
    `FIXED_VALUE_BOUNDARY`.

    Parameters
    ----------
    shape : tuple
        Shape of the structured grid
    node_status : array_like, optional
        Status of each node in the grid.
    return_count : boolean, optional
        If `True`, also return an array of active link counts per node.

    Returns
    -------
    links : (2, N) ndarray
        Active link IDs for each node.
    count : ndarray
        Number of active links per node.

    Examples
    --------
    Get the active link IDs for a grid of 3 nodes by 4 nodes. The first row
    lists links entering nodes from the bottom, and the second links entering
    from the left.

    >>> from landlab.utils.structured_grid import setup_active_inlink_matrix2
    >>> setup_active_inlink_matrix2((3, 4), return_count=False)
    array([[-1, -1, -1, -1, -1,  4,  5, -1, -1, 11, 12, -1],
           [-1, -1, -1, -1, -1,  7,  8,  9, -1, -1, -1, -1]])
    >>> _, count = setup_active_inlink_matrix2((3, 4))
    >>> count
    array([0, 0, 0, 0, 0, 2, 2, 1, 0, 1, 1, 0])
    """
    links = active_inlinks2(shape, node_status=node_status)
    if return_count:
        return links, active_inlink_count_per_node(shape)
    else:
        return links


def node_index_with_halo(shape, halo_indices=BAD_INDEX_VALUE):
    """Array of links with a halo of no-data values.

    Examples
    --------
    >>> from landlab.utils.structured_grid import node_index_with_halo
    >>> node_index_with_halo((2, 3), halo_indices=-1)
    array([[-1, -1, -1, -1, -1],
           [-1,  0,  1,  2, -1],
           [-1,  3,  4,  5, -1],
           [-1, -1, -1, -1, -1]])
    """
    shape_with_halo = np.array(shape) + 2

    ids = np.empty(shape_with_halo, dtype=np.int)

    (interiors, boundaries) = (
        interior_nodes(shape_with_halo),
        boundary_nodes(shape_with_halo),
    )

    ids.flat[interiors] = range(interior_node_count(shape_with_halo))
    ids.flat[boundaries] = halo_indices

    return ids


def cell_index_with_halo(shape, halo_indices=BAD_INDEX_VALUE, inactive_indices=None):
    """Array of cells with a halo of no-data values.

    Examples
    --------
    >>> from landlab.utils.structured_grid import cell_index_with_halo
    >>> cell_index_with_halo((2, 3), halo_indices=-1)
    array([[-1, -1, -1, -1, -1],
           [-1,  0,  1,  2, -1],
           [-1,  3,  4,  5, -1],
           [-1, -1, -1, -1, -1]])

    >>> cell_index_with_halo((2, 3), halo_indices=-1, inactive_indices=-1)
    array([[-1, -1, -1, -1, -1],
           [-1, -1, -1, -1, -1],
           [-1, -1, -1, -1, -1],
           [-1, -1, -1, -1, -1]])
    """
    ids = node_index_with_halo(shape, halo_indices=halo_indices)
    if inactive_indices is not None:
        ids[:, (1, -2)] = inactive_indices
        ids[(1, -2), :] = inactive_indices

    return ids


def _neighbor_node_ids(ids_with_halo):
    """Matrix of four neighbor nodes for each node."""
    shape = (ids_with_halo.shape[0] - 2, ids_with_halo.shape[1] - 2)
    kwds = {
        "strides": ids_with_halo.strides,
        "buffer": ids_with_halo,
        "dtype": ids_with_halo.dtype,
    }

    kwds["offset"] = ids_with_halo.itemsize * (ids_with_halo.shape[1])
    west_ids = np.ndarray(shape, **kwds)

    kwds["offset"] = ids_with_halo.itemsize * (ids_with_halo.shape[1] + 2)
    east_ids = np.ndarray(shape, **kwds)

    kwds["offset"] = ids_with_halo.itemsize
    south_ids = np.ndarray(shape, **kwds)

    kwds["offset"] = ids_with_halo.itemsize * (ids_with_halo.shape[1] * 2 + 1)
    north_ids = np.ndarray(shape, **kwds)

    return np.vstack((east_ids.flat, north_ids.flat, west_ids.flat, south_ids.flat))


def _centered_node_ids(ids_with_halo):
    """Array of nodes taken from a matrix of nodes with a halo."""
    shape = (ids_with_halo.shape[0] - 2, ids_with_halo.shape[1] - 2)
    kwds = {
        "strides": ids_with_halo.strides,
        "buffer": ids_with_halo,
        "dtype": ids_with_halo.dtype,
    }

    kwds["offset"] = ids_with_halo.itemsize * (ids_with_halo.shape[1] + 1)
    return np.ndarray(shape, **kwds)


def neighbor_node_ids(shape, inactive=BAD_INDEX_VALUE):
    """Matrix of four neighbor nodes for each node."""
    return linked_neighbor_node_ids(shape, [], inactive=inactive)


def linked_neighbor_node_ids(
    shape, closed_boundary_nodes, open_boundary_nodes=None, inactive=BAD_INDEX_VALUE
):
    """Matrix of four neighbor nodes for each node."""
    if open_boundary_nodes is None:
        open_boundary_nodes = []

    ids_with_halo = node_index_with_halo(shape, halo_indices=inactive)

    # Everything that touches a closed boundary is inactive
    if len(closed_boundary_nodes) > 0:
        ids = _centered_node_ids(ids_with_halo)
        ids.flat[closed_boundary_nodes] = inactive

    neighbors = _neighbor_node_ids(ids_with_halo)

    # Everything that a closed boundary touches is inactive
    if len(closed_boundary_nodes) > 0:
        neighbors[:, closed_boundary_nodes] = inactive

    if len(open_boundary_nodes) > 0:
        _set_open_boundary_neighbors(neighbors, open_boundary_nodes, inactive)

    return neighbors


def _set_open_boundary_neighbors(neighbors, open_boundary_nodes, value):
    """Set values for open-boundary neighbor-nodes."""
    open_boundary_neighbors = neighbors[:, open_boundary_nodes]
    is_open_boundary_neighbor = _find_open_boundary_neighbors(
        neighbors, open_boundary_nodes
    )
    nodes = np.choose(is_open_boundary_neighbor, (open_boundary_neighbors, value))
    neighbors[:, open_boundary_nodes] = nodes


def _find_open_boundary_neighbors(neighbors, open_boundary_nodes):
    """Array of booleans that indicate if a neighbor is an open boundary."""
    open_boundary_neighbors = neighbors[:, open_boundary_nodes]
    is_open_boundary_neighbor = np.in1d(open_boundary_neighbors, open_boundary_nodes)
    is_open_boundary_neighbor.shape = (neighbors.shape[0], len(open_boundary_nodes))
    return is_open_boundary_neighbor


def neighbor_node_array(shape, **kwds):
    """Array of neighbor nodes.

    Examples
    --------
    >>> from landlab.utils.structured_grid import neighbor_node_array
    >>> neighbors = neighbor_node_array((2, 3), inactive=-1)
    >>> neighbors.T
    array([[ 1,  3, -1, -1],
           [ 2,  4,  0, -1],
           [-1,  5,  1, -1],
           [ 4, -1, -1,  0],
           [ 5, -1,  3,  1],
           [-1, -1,  4,  2]])
    """
    closed_boundary_nodes = kwds.pop("closed_boundary_nodes", [])
    open_boundary_nodes = kwds.get("open_boundary_nodes", [])

    if len(closed_boundary_nodes) > 0 or len(open_boundary_nodes):
        neighbors = linked_neighbor_node_ids(shape, closed_boundary_nodes, **kwds)
    else:
        neighbors = neighbor_node_ids(shape, **kwds)

    return neighbors


def neighbor_cell_array(shape, out_of_bounds=BAD_INDEX_VALUE, contiguous=True):
    """Array of neighbor cells.

    Examples
    --------
    >>> from landlab.utils.structured_grid import neighbor_cell_array
    >>> neighbors = neighbor_cell_array((2, 3), out_of_bounds=-1)
    >>> len(neighbors) == 0
    True

    >>> neighbors = neighbor_cell_array((3, 3), out_of_bounds=-1)
    >>> neighbors
    array([[-1, -1, -1, -1]])

    >>> neighbors = neighbor_cell_array((5, 4), out_of_bounds=-1)
    >>> neighbors # doctest: +NORMALIZE_WHITESPACE
    array([[ 1,  2, -1, -1], [-1,  3,  0, -1],
           [ 3,  4, -1,  0], [-1,  5,  2,  1],
           [ 5, -1, -1,  2], [-1, -1,  4,  3]])
    """
    if cell_count(shape) > 0:
        shape = np.array(shape) - 2
        ids = node_index_with_halo(shape, halo_indices=out_of_bounds)

        neighbors = np.vstack(
            (
                ids[1 : shape[0] + 1, 2:].flat,
                ids[2:, 1 : shape[1] + 1].flat,
                ids[1 : shape[0] + 1, : shape[1]].flat,
                ids[: shape[0], 1 : shape[1] + 1].flat,
            )
        ).T
        if contiguous:
            return neighbors.copy()
        else:
            return neighbors
    else:
        return np.array([], dtype=np.int)


def diagonal_node_array(
    shape, out_of_bounds=BAD_INDEX_VALUE, contiguous=True, boundary_node_mask=None
):
    """Array of diagonal nodes.

    Creates a list of IDs of the diagonal cells to each cell, as a 2D array.
    Only interior cells are assigned neighbors; boundary cells get -1 for
    each neighbor.  The order of the diagonal cells is [topright, topleft,
    bottomleft, bottomright].

    Examples
    --------
    >>> from landlab.utils.structured_grid import diagonal_node_array
    >>> diags = diagonal_node_array((2, 3), out_of_bounds=-1)
    >>> diags
    array([[ 4, -1, -1, -1],
           [ 5,  3, -1, -1],
           [-1,  4, -1, -1],
           [-1, -1, -1,  1],
           [-1, -1,  0,  2],
           [-1, -1,  1, -1]])
    >>> diags.flags['C_CONTIGUOUS']
    True
    >>> diags = diagonal_node_array((2, 3), out_of_bounds=-1, contiguous=False)
    >>> diags.flags['C_CONTIGUOUS']
    False
    """
    # NG didn't touch this, but she thinks this should be nodes, not cells.
    ids = node_index_with_halo(shape, halo_indices=out_of_bounds)

    diags = np.vstack(
        (
            ids[2:, 2:].flat,
            ids[2:, : shape[1]].flat,
            ids[: shape[0], : shape[1]].flat,
            ids[: shape[0], 2:].flat,
        )
    ).T

    if boundary_node_mask is not None:
        boundaries = np.empty(4, dtype=np.int)
        boundaries.fill(boundary_node_mask)
        diags[boundary_nodes(shape)] = boundaries

    if contiguous:
        return diags.copy()
    else:
        return diags


def diagonal_cell_array(shape, out_of_bounds=BAD_INDEX_VALUE, contiguous=True):
    """Array of diagonal cells.

    Construct a matrix of cell indices to each diagonally adjacent cell of a
    structured grid. If a cell does not have a diagonal neighbor, set the
    index for that neighbor to *out_of_bounds*.

    Examples
    --------
    A grid without any cells returns an empty array.

    >>> from landlab.utils.structured_grid import diagonal_cell_array
    >>> diags = diagonal_cell_array((2, 3), out_of_bounds=-1)
    >>> len(diags) == 0
    True

    A grid that has only one cell does not have any neighbors so all of its
    diagonals are set to *out_of_bounds*.
    >>> diags = diagonal_cell_array((3, 3), out_of_bounds=-1)
    >>> diags
    array([[-1, -1, -1, -1]])

    >>> diags = diagonal_cell_array((4, 4), out_of_bounds=-1)
    >>> diags # doctest: +NORMALIZE_WHITESPACE
    array([[ 3, -1, -1, -1], [-1,  2, -1, -1],
           [-1, -1, -1,  1], [-1, -1,  0, -1]])

    >>> diags = diagonal_cell_array((4, 5), out_of_bounds=-1)
    >>> diags # doctest: +NORMALIZE_WHITESPACE
    array([[ 4, -1, -1, -1], [ 5,  3, -1, -1], [-1,  4, -1, -1],
           [-1, -1, -1,  1], [-1, -1,  0,  2], [-1, -1,  1, -1]])
    """
    if cell_count(shape) > 0:
        shape = np.array(shape) - 2
        ids = node_index_with_halo(shape, halo_indices=out_of_bounds)

        diags = np.vstack(
            (
                ids[2:, 2:].flat,
                ids[2:, : shape[1]].flat,
                ids[: shape[0], : shape[1]].flat,
                ids[: shape[0], 2:].flat,
            )
        ).T
        if contiguous:
            return diags.copy()
        else:
            return diags
    else:
        return np.array([], dtype=np.int)


def node_has_boundary_neighbor(neighbors, diagonals, out_of_bounds=BAD_INDEX_VALUE):
    """Array of booleans that indicate if a node has a boundary neighbor.

    .. note::

        DEJH thinks this method is broken since terminology update: it returns
        closed neighbors, not boundary neighbors.
    """
    return out_of_bounds in neighbors | out_of_bounds in diagonals


def reshape_array(shape, array, flip_vertically=False, copy=False):
    """Reshape a flat array.

    Examples
    --------
    >>> from landlab.utils.structured_grid import reshape_array
    >>> x = np.arange(12.)
    >>> y = reshape_array((3, 4), x)
    >>> y.shape == (3, 4)
    True
    >>> y
    array([[  0.,   1.,   2.,   3.],
           [  4.,   5.,   6.,   7.],
           [  8.,   9.,  10.,  11.]])
    >>> y.flags['C_CONTIGUOUS']
    True
    >>> x[0] = -1
    >>> y[0, 0]
    -1.0

    >>> x = np.arange(12.)
    >>> y = reshape_array((3, 4), x, flip_vertically=True)
    >>> y
    array([[  8.,   9.,  10.,  11.],
           [  4.,   5.,   6.,   7.],
           [  0.,   1.,   2.,   3.]])
    >>> y.flags['C_CONTIGUOUS']
    False
    >>> x[0] = -1
    >>> y[-1, 0]
    -1.0
    """
    reshaped_array = array.view()

    try:
        reshaped_array.shape = shape
    except ValueError:
        raise

    if flip_vertically:
        flipped_array = reshaped_array[::-1, :]
        if copy:
            return flipped_array.copy()
        else:
            return flipped_array
    else:
        if copy:
            return reshaped_array.copy()
        else:
            return reshaped_array


def nodes_around_points_on_unit_grid(shape, coords, mode="raise"):
    """Array of nodes around x, y points on a grid of unit spacing.

    Returns the nodes around a point on a structured grid with unit spacing
    and zero origin.

    Examples
    --------
    >>> from landlab.utils.structured_grid import (
    ...     nodes_around_points_on_unit_grid)
    >>> nodes_around_points_on_unit_grid((3, 3), (.1, .1))
    array([0, 3, 4, 1])

    >>> nodes_around_points_on_unit_grid((3, 3), (1., 1.))
    array([4, 7, 8, 5])
    """
    if isinstance(coords[0], np.ndarray):
        (rows, cols) = (as_id_array(coords[0]), as_id_array(coords[1]))
    else:
        (rows, cols) = (int(coords[0]), int(coords[1]))

    return as_id_array(
        np.ravel_multi_index(
            ((rows, rows + 1, rows + 1, rows), (cols, cols, cols + 1, cols + 1)),
            shape,
            mode=mode,
        ).T
    )


def nodes_around_points(shape, coords, spacing=(1.0, 1.0), origin=(0.0, 0.0)):
    """Array of nodes around x, y points on a grid of non-unit spacing.

    Returns the nodes around a point on a structured grid with row and column
    *spacing*, and *origin*.

    Examples
    --------
    >>> from landlab.utils.structured_grid import nodes_around_points
    >>> x = np.array([.9, 1.])
    >>> y = np.array([.1, 1.])
    >>> nodes_around_points((3, 3), (y, x))
    array([[0, 3, 4, 1],
           [4, 7, 8, 5]])

    >>> nodes_around_points((3, 3), (2., 1.))
    Traceback (most recent call last):
        ...
    ValueError: invalid entry in coordinates array
    """
    return as_id_array(
        nodes_around_points_on_unit_grid(
            shape,
            (
                (coords[0] - origin[0]) / spacing[0],
                (coords[1] - origin[1]) / spacing[1],
            ),
        )
    )


def nodes_around_point(shape, coords, spacing=(1.0, 1.0)):
    """Array of nodes around a single point on a grid of non-unit spacing."""
    node_id = int(coords[0] // spacing[0] * shape[1] + coords[1] // spacing[1])
    if node_id + shape[1] + 1 >= shape[0] * shape[1] or node_id < 0:
        raise ValueError("invalid entry in coordinates array")

    return np.array([node_id, node_id + shape[1], node_id + shape[1] + 1, node_id + 1])

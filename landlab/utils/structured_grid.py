#! /usr/bin/env python
import numpy as np
import itertools

from ..grid.base import (CORE_NODE, FIXED_GRADIENT_BOUNDARY,
                         FIXED_VALUE_BOUNDARY, TRACKS_CELL_BOUNDARY,
                         CLOSED_BOUNDARY, BAD_INDEX_VALUE)


def node_count(shape):
    """
    The total number of nodes in a structured grid with dimensions given
    by the tuple, *shape*. Where *shape* is the number of node rows and
    node columns.

    >>> from landlab.utils.structured_grid import node_count
    >>> node_count((3, 4))
    12
    """
    assert(len(shape) == 2)
    return shape[0] * shape[1]


def interior_node_count(shape):
    """
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
    assert(len(shape) == 2)
    try:
        assert(np.min(shape) > 2)
    except AssertionError:
        return 0
    else:
        return (shape[0] - 2) * (shape[1] - 2)


def cell_count(shape):
    """
    The total number of cells in a structured grid with dimensions, *shape*.
    Where *shape* is a tuple that gives the dimensions of the grid as number
    of rows of nodes followed by number of columns of nodes.

    >>> from landlab.utils.structured_grid import cell_count
    >>> cell_count((3, 4))
    2
    >>> cell_count((1, 4))
    0
    """
    assert(len(shape) == 2)

    if np.min(shape) > 2:
        return (shape[0] - 2) * (shape[1] - 2)
    else:
        return 0


def active_cell_count(shape):
    """
    Number of active cells. By default, all cells are active so this is
    the same as cell_count. (active = core+open boundary)
    """
    return cell_count(shape)

def core_cell_count(shape):
    """
    Number of core cells. By default, all cells are core so this is
    the same as cell_count.
    """
    return cell_count(shape)

def active_link_count(shape):
    """
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
    assert(len(shape) == 2)
    try:
        assert(np.min(shape) > 2)
    except AssertionError:
        return 0
    else:
        return 2 * shape[0] * shape[1] - 3 * (shape[0] + shape[1]) + 4


def link_count(shape):
    """
    Total (active and inactive) number of links in a structured grid with
    dimensions, *shape*. This is the number of to-links and from-links, not
    the total of the two.

    >>> from landlab.utils.structured_grid import link_count
    >>> link_count((3,2))
    7
    """
    assert(len(shape) == 2)
    return shape[1] * (shape[0] - 1) + shape[0] * (shape[1] - 1)


def vertical_link_count(shape):
    assert(len(shape) == 2)
    return (shape[0] - 1) * shape[1]


def horizontal_link_count(shape):
    assert(len(shape) == 2)
    return shape[0] * (shape[1] - 1)


def boundary_cell_count(shape):
    """
    .. deprecated:: 0.6
    Deprecated due to imprecise terminology; this function makes little sense.
    Use perimeter_node_count() instead.
    Number of cells that are on the perimeter of a structured grid with
    dimensions, *shape*, and thus boundary cells. In fact, cells centered on
    boundary nodes are not really cells. If they were, though, this is how
    many there would be.

    **** Shouldn't be deprecated. This routine returns the cells on the
    boundary. Not the cells surrounding boundary nodes because there aren't
    cells around boundary nodes by definition as previously understood.
      - SN  30Nov14  ****

    Examples
    --------
    >>> from landlab.utils.structured_grid import boundary_cell_count
    >>> boundary_cell_count((3, 4))
    10
    """
    assert(len(shape) == 2)
    return 2 * (shape[0] - 2) + 2 * (shape[1] - 2) + 4


def perimeter_node_count(shape):
    """
    Number of nodes that are on the perimeter of a structured grid with
    dimensions, *shape*, and thus boundary nodes.

    >>> from landlab.utils.structured_grid import perimeter_node_count
    >>> perimeter_node_count((3, 4))
    10
    """
    assert(len(shape) == 2)
    return 2 * (shape[0] - 2) + 2 * (shape[1] - 2) + 4


def interior_cell_count(shape):
    """
    Number of interior cells. Since cells are only defined on interior nodes,
    this is the same as cell_count.
    """
    return cell_count(shape)


def face_count(shape):
    """
    Total number of faces in a structured grid with dimensions, *shape*. Each
    cell has four faces, and shared faces only count once.

    >>> from landlab.utils.structured_grid import face_count
    >>> face_count((3, 4))
    7
    """
    assert(len(shape) == 2)
    if np.min(shape) > 2:
        return ((shape[0] - 1) * (shape[1] - 2) +
                (shape[0] - 2) * (shape[1] - 1))
    else:
        return 0


def active_face_count(shape):
    """
    Total number of active faces in a structured grid with dimensions,
    *shape*. Each cell has four faces, and shared faces only count once.
    An active face is one that has a corresponing active link.

    >>> from landlab.utils.structured_grid import active_face_count
    >>> active_face_count((3, 4))
    7
    """
    return face_count(shape)


def top_index_iter(shape):
    """
    Iterator for the top boundary indices of a structured grid.
    """
    return xrange(shape[1] * (shape[0] - 1), shape[0] * shape[1])


def bottom_index_iter(shape):
    """
    Iterator for the bottom boundary indices of a structured grid.
    """
    return xrange(0, shape[1])


def left_index_iter(shape):
    """
    Iterator for the left boundary indices of a structured grid.
    """
    return xrange(0, shape[0] * shape[1], shape[1])


def right_index_iter(shape):
    """
    Iterator for the right boundary indices of a structured grid.
    """
    return xrange(shape[1] - 1, shape[0] * shape[1], shape[1])


def left_right_iter(shape, *args):
    """
    Iterator for the left and right boundary indices of a structured grid.
    This iterates over the indices in order rather than iterating all of
    the left boundary and then all of the right boundary.

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
        iter_rows = xrange(0, shape[0], 1)
    elif len(args) == 1:
        iter_rows = xrange(0, args[0], 1)
    elif len(args) == 2:
        iter_rows = xrange(args[0], args[1], 1)
    elif len(args) == 3:
        iter_rows = xrange(args[0], args[1], args[2])

    for row in iter_rows:
        yield row * shape[1]
        yield row * shape[1] + shape[1] - 1


def bottom_top_iter(shape):
    """
    Iterates of the bottom indices and then the top indices of a structured
    grid.

    >>> import numpy as np
    >>> from landlab.utils.structured_grid import bottom_top_iter
    >>> np.fromiter(bottom_top_iter((4, 3)), dtype=np.int)
    array([ 0,  1,  2,  9, 10, 11])
    """
    return itertools.chain(bottom_index_iter(shape),
                           top_index_iter(shape))


def boundary_iter(shape):
    """
    .. deprecated:: 0.6
    Deprecated due to imprecise terminology. This is really perimeter_iter
    (see below).
    Iterates over all of the boundary node indices of a structured grid in
    order.

    >>> import numpy as np
    >>> from landlab.utils.structured_grid import boundary_iter
    >>> np.fromiter(boundary_iter((4, 3)), dtype=np.int)
    array([ 0,  1,  2,  3,  5,  6,  8,  9, 10, 11])
    """
    return itertools.chain(bottom_index_iter(shape),
                           left_right_iter(shape, 1, shape[0] - 1),
                           top_index_iter(shape))

def perimeter_iter(shape):
    """
    Iterates over all of the perimeter node indices of a structured grid in
    order.

    >>> import numpy as np
    >>> from landlab.utils.structured_grid import perimeter_iter
    >>> np.fromiter(perimeter_iter((4, 3)), dtype=np.int)
    array([ 0,  1,  2,  3,  5,  6,  8,  9, 10, 11])
    """
    return itertools.chain(bottom_index_iter(shape),
                           left_right_iter(shape, 1, shape[0] - 1),
                           top_index_iter(shape))

def boundary_nodes(shape):
    """
    .. deprecated:: 0.6
    Deprecated due to imprecise terminology. This is really perimeter_iter
    (see below).
    An array of the indices of the boundary nodes.

    >>> from landlab.utils.structured_grid import boundary_nodes
    >>> boundary_nodes((3, 4))
    array([ 0,  1,  2,  3,  4,  7,  8,  9, 10, 11])
    """
    return np.fromiter(boundary_iter(shape), dtype=np.int)

def perimeter_nodes(shape):
    """
    An array of the indices of the perimeter nodes of a structured grid.

    >>> from landlab.utils.structured_grid import perimeter_nodes
    >>> perimeter_nodes((3, 4))
    array([ 0,  1,  2,  3,  4,  7,  8,  9, 10, 11])
    """
    return np.fromiter(perimeter_iter(shape), dtype=np.int)

def corners(shape):
    """
    An array of the indices of the grid corner nodes.
    """
    return np.array([0, shape[1]-1, shape[1]*(shape[0]-1), shape[1]*shape[0]-1])

def bottom_edge_node_ids(shape):
    return np.fromiter(bottom_index_iter(shape), dtype=np.int)


def top_edge_node_ids(shape):
    return np.fromiter(top_index_iter(shape), dtype=np.int)


def left_edge_node_ids(shape):
    return np.fromiter(left_index_iter(shape), dtype=np.int)


def right_edge_node_ids(shape):
    return np.fromiter(right_index_iter(shape), dtype=np.int)


def interior_iter(shape):
    """
    Iterate over the interior nodes of a structured grid.

    >>> import numpy as np
    >>> from landlab.utils.structured_grid import interior_iter
    >>> np.fromiter(interior_iter((4, 3)), dtype=np.int)
    array([4, 7])
    """
    interiors = []
    interiors_per_row = shape[1] - 2
    for row in xrange(shape[1] + 1, shape[1] * (shape[0] - 1), shape[1]):
        interiors.append(xrange(row , row + interiors_per_row))
    return itertools.chain(*interiors)


def interior_nodes(shape):
    return np.fromiter(interior_iter(shape), dtype=np.int)


def node_coords(shape, *args):
    """
    Get x, y coordinates for nodes in a structured grid with dimensions,
    *shape*. Use the optional argument *spacing* to give the spacing in each
    dimension, and *origin* the start of the coordinates in each dimension.

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
        assert(len(spacing) == len(shape))

    try:
        origin = args[1]
    except IndexError:
        origin = np.zeros(len(shape), dtype=np.float)
    else:
        assert(len(origin) == len(origin))

    node_count = np.prod(shape)

    row_y = np.arange(shape[0]) * spacing[0] + origin[0]
    col_x = np.arange(shape[1]) * spacing[1] + origin[1]

    (node_x, node_y) = np.meshgrid(col_x, row_y)

    node_x.shape = (node_count, )
    node_y.shape = (node_count, )

    return (node_x, node_y)


def active_cell_index(shape):
    """
    For many instances, core_cell_index() may be preferred.
    Ordered indices of the active (core+open boundary) cells of a structured
    grid.
    """
    return np.arange(active_cell_count(shape))


def core_cell_index(shape):
    """
    Ordered indices of the core cells of a structured grid.
    """
    return np.arange(core_cell_count(shape))


def active_cell_node(shape):
    """
    For many instances, core_cell_node() may be preferred.
    Indices of the nodes belonging to each active (core + open boundary) cell.
    Since all cells are active in the default case, this is the same as
    node_index_at_cells.

    >>> from landlab.utils.structured_grid import active_cell_node
    >>> active_cell_node((4,3))
    array([4, 7])
    """
    return node_index_at_cells(shape)


def core_cell_node(shape):
    """
    Indices of the nodes belonging to each core cell.
    Since all cells are core in the default case, this is the same as
    node_index_at_cells.

    >>> from landlab.utils.structured_grid import core_cell_node
    >>> core_cell_node((4,3))
    array([4, 7])
    """
    return node_index_at_cells(shape)


def active_cell_index_at_nodes(shape, boundary_node_index=BAD_INDEX_VALUE):
    """
    .. deprecated:: 0.6
    Up to date, unambiguous terminology core_cell_index_at_nodes() preferred.
    Indices of the active cells associated with the nodes of the structured grid.
    For nodes that don't have a cell (that is, boundary nodes) set indices
    to BAD_INDEX_VALUE. Use the *boundary_node_index* keyword to change
    the value of indices to boundary nodes.

    Note that all three functions [X_]cell_index_at_nodes are equivalent.

    >>> from landlab.utils.structured_grid import active_cell_index_at_nodes
    >>> active_cell_index_at_nodes((3, 4), boundary_node_index=-1) # doctest: +NORMALIZE_WHITESPACE
    array([-1, -1, -1, -1,
           -1,  0,  1, -1,
           -1, -1, -1, -1])
    """
    n_nodes = node_count(shape)

    node_ids = np.empty(n_nodes, dtype=np.int)

    node_ids[boundary_nodes(shape)] = boundary_node_index
    node_ids[interior_nodes(shape)] = np.arange(interior_node_count(shape))

    return node_ids

def core_cell_index_at_nodes(shape, boundary_node_index=BAD_INDEX_VALUE):
    """
    Indices of the core cells associated with the nodes of the structured grid.
    For nodes that don't have a cell (that is, boundary nodes) set indices
    to BAD_INDEX_VALUE. Use the *boundary_node_index* keyword to change
    the value of indices to boundary nodes.

    Note that all three functions [X_]cell_index_at_nodes are equivalent.

    >>> from landlab.utils.structured_grid import core_cell_index_at_nodes
    >>> core_cell_index_at_nodes((3, 4), boundary_node_index=-1) # doctest: +NORMALIZE_WHITESPACE
    array([-1, -1, -1, -1,
           -1,  0,  1, -1,
           -1, -1, -1, -1])
    """
    n_nodes = node_count(shape)

    node_ids = np.empty(n_nodes, dtype=np.int)

    node_ids[perimeter_nodes(shape)] = boundary_node_index
    node_ids[interior_nodes(shape)] = np.arange(interior_node_count(shape))

    return node_ids

def cell_index_at_nodes(shape, boundary_node_index=BAD_INDEX_VALUE):
    """
    Indices of the core cells associated with the nodes of the structured grid.
    For nodes that don't have a cell (that is, boundary nodes) set indices
    to BAD_INDEX_VALUE. Use the *boundary_node_index* keyword to change
    the value of indices to boundary nodes.

    Note that all three functions [X_]cell_index_at_nodes are equivalent.

    >>> from landlab.utils.structured_grid import cell_index_at_nodes
    >>> cell_index_at_nodes((3, 4), boundary_node_index=-1) # doctest: +NORMALIZE_WHITESPACE
    array([-1, -1, -1, -1,
           -1,  0,  1, -1,
           -1, -1, -1, -1])
    """
    n_nodes = node_count(shape)

    node_ids = np.empty(n_nodes, dtype=np.int)

    node_ids[perimeter_nodes(shape)] = boundary_node_index
    node_ids[interior_nodes(shape)] = np.arange(interior_node_count(shape))

    return node_ids


def node_index_at_cells(shape):
    """
    Indices of the nodes belonging to each cell.

    >>> from landlab.utils.structured_grid import node_index_at_cells
    >>> node_index_at_cells((4, 3))
    array([4, 7])
    """
    node_ids = np.arange(node_count(shape))
    node_ids.shape = shape

    cell_node = node_ids[1:-1, 1:-1].copy()
    cell_node.shape = ((shape[0] - 2) * (shape[1] - 2), )

    return cell_node


def node_index_at_link_ends(shape):
    node_ids = np.arange(np.prod(shape))
    node_ids.shape = shape

    return (node_index_at_link_head(node_ids),
            node_index_at_link_tail(node_ids))


def inlink_index_at_node(shape):
    return inlinks(shape, return_count=False)


def outlink_index_at_node(shape):
    return outlinks(shape, return_count=False)


def node_index_at_link_tail(node_ids):
    vertical_links = node_ids[1:, :]
    horizontal_links = node_ids[:, 1:]
    return np.concatenate((vertical_links.flat, horizontal_links.flat))


def node_index_at_link_head(node_ids):
    vertical_links = node_ids[:-1, :]
    horizontal_links = node_ids[:, :-1]
    return np.concatenate((vertical_links.flat, horizontal_links.flat))


def face_index_at_links(shape, actives=None,
                        inactive_link_index=BAD_INDEX_VALUE):
    """
    Returns an array that maps link ids to face ids. For inactive links,
    which do not have associated faces, set their ids to
    *inactive_link_index*. Use the *actives* keyword to specify an array that
    contains the ids of all active links in the grid. The default assumes
    that only the perimeter nodes are inactive.


    >>> from landlab.utils.structured_grid import face_index_at_links
    >>> faces = face_index_at_links((3, 4), inactive_link_index=-1)
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


def node_status(shape, boundary_status=FIXED_VALUE_BOUNDARY):
    """
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

    Return the link IDs for links that are *active* in a structured grid of quadrilaterals. Use the
    *node_status_array* keyword to specify the status for each of the grid's nodes. If not given, each of the
    perimeter nodes is assumed to be `FIXED_VALUE_BOUNDARY`.

    Use the *link_nodes* keyword to provide, as a tuple of arrays, that give the *from-node* and the *to-node*
    for each for each link in the grid.

    Parameters
    ----------
    shape : tuple
        Shape of grid as number of node rows and columns.
    node_status_array : array_like, optional
        Status of each grid node.
    link_nodes : array_like, optional

    Examples
    --------
    Because, by default, the perimeter nodes are `FIXED_VALUE_BOUNDARY` nodes, only links attached to the
    interior nodes are *active*.

    >>> from landlab.utils.structured_grid import active_links
    >>> from landlab import CLOSED_BOUNDARY, CORE_NODE
    >>> active_links((3, 4))
    array([ 1,  2,  5,  6, 11, 12, 13])

    If all the perimeter nodes `CLOSED_BOUNDARY` nodes, the only active link is between the two core nodes.

    >>> node_status = np.ones(3 * 4) * CLOSED_BOUNDARY
    >>> node_status[5:7] = CORE_NODE
    >>> active_links((3, 4), node_status_array=node_status)
    array([12])

    You can also provide a list of all the *from_nodes* and *to_nodes* for the grid. The following describes
    a grid with only a single link (between nodes 5 and 6).

    >>> active_links((3, 4), link_nodes=(np.array([5]), np.array([6])))
    array([0])
    """
    if node_status_array is None:
        node_status_array = node_status(shape)

    if link_nodes is None:
        (link_from_node, link_to_node) = node_index_at_link_ends(shape)
    else:
        (link_from_node, link_to_node) = link_nodes

    from_node_status = node_status_array[link_from_node]
    to_node_status = node_status_array[link_to_node]

    active_links = (((from_node_status == CORE_NODE) & ~
                     (to_node_status == CLOSED_BOUNDARY)) |
                    ((to_node_status == CORE_NODE) & ~
                     (from_node_status == CLOSED_BOUNDARY)))

    (active_links, ) = np.where(active_links)

    return active_links


def active_face_index(shape):
    return np.arange(active_face_count(shape))


def inlinks(shape, tonodes=None, return_count=True):
    links = np.vstack((south_links(shape), west_links(shape)))
    links.shape = (2, node_count(shape))
    return links


def outlinks(shape, tonodes=None, return_count=True):
    links = np.vstack((north_links(shape), east_links(shape)))
    links.shape = (2, node_count(shape))
    return links


def active_inlinks(shape, node_status=None, return_count=True):
    links = np.vstack((active_south_links(shape, node_status),
                       active_west_links(shape, node_status)))
    links.shape = (2, node_count(shape))
    return links


def active_outlinks(shape, node_status=None, return_count=True):
    links = np.vstack((active_north_links(shape, node_status),
                       active_east_links(shape, node_status)))
    links.shape = (2, node_count(shape))
    return links


def vertical_link_ids(shape):
    link_ids = np.arange(0, vertical_link_count(shape), dtype=np.int)
    link_ids.shape = (shape[0] - 1, shape[1])
    return link_ids


def horizontal_link_ids(shape):
    link_id_offset = vertical_link_count(shape)
    link_ids = np.arange(link_id_offset,
                         link_id_offset + horizontal_link_count(shape),
                         dtype=np.int)
    link_ids.shape = (shape[0], shape[1] - 1)
    return link_ids


def vertical_active_link_count(shape, node_status=None):
    if node_status is not None:
        is_inactive = (node_status == 0)
        is_inactive.shape = shape

        inactive_outlinks = is_inactive[:-1, 1:-1]
        inactive_inlinks = is_inactive[1:, 1:-1]
        inactive_links = inactive_outlinks | inactive_inlinks

        return (shape[0] - 1) * (shape[1] - 2) - np.sum(inactive_links.flat)
    else:
        return (shape[0] - 1) * (shape[1] - 2)


def horizontal_active_link_count(shape, node_status=None):
    if node_status is not None:
        is_inactive = (node_status == 0)
        is_inactive.shape = shape

        inactive_outlinks = is_inactive[1:-1, :-1]
        inactive_inlinks = is_inactive[1:-1, 1:]
        inactive_links = inactive_outlinks | inactive_inlinks

        return (shape[0] - 2) * (shape[1] - 1) - np.sum(inactive_links.flat)
    else:
        return (shape[0] - 2) * (shape[1] - 1)


def vertical_inactive_link_mask(shape, node_status):
    is_inactive = (node_status == 0)
    is_inactive.shape = shape

    inactive_outlinks = is_inactive[:-1, 1:-1]
    inactive_inlinks = is_inactive[1:, 1:-1]
    return inactive_outlinks | inactive_inlinks


def horizontal_inactive_link_mask(shape, node_status):
    is_inactive = (node_status == 0)
    is_inactive.shape = shape

    inactive_outlinks = is_inactive[1:-1, :-1]
    inactive_inlinks = is_inactive[1:-1, 1:]
    return inactive_outlinks | inactive_inlinks


#def vertical_active_link_ids(shape):
#    link_ids = np.arange(vertical_active_link_count(shape), dtype=np.int)
#    link_ids.shape = (shape[0] - 1, shape[1] - 2)
#    return link_ids


def vertical_active_link_ids(shape, node_status=None):
    if node_status is None:
        link_ids = np.arange(vertical_active_link_count(shape), dtype=np.int)
        #link_ids.shape = (shape[0] - 1, shape[1] - 2)
        #return link_ids
    else:
        inactive_links = vertical_inactive_link_mask(shape, node_status)
        inactive_links.shape = (inactive_links.size, )
        active_link_count = inactive_links.size - np.sum(inactive_links)

        link_ids = np.empty(inactive_links.size)
        link_ids[inactive_links] = - 1
        link_ids[~ inactive_links] = np.arange(active_link_count, dtype=int)

    link_ids.shape = (shape[0] - 1, shape[1] - 2)
    return link_ids


def horizontal_active_link_ids(shape, node_status=None):
    if node_status is None:
        link_id_offset = vertical_active_link_count(shape)
        link_ids = np.arange(
            link_id_offset,
            link_id_offset + horizontal_active_link_count(shape),
            dtype=np.int)
    else:
        link_id_offset = vertical_active_link_count(shape,
                                                    node_status=node_status)
        inactive_links = horizontal_inactive_link_mask(shape, node_status)
        inactive_links.shape = (inactive_links.size, )
        active_link_count = inactive_links.size - np.sum(inactive_links)

        link_ids = np.empty(inactive_links.size)
        link_ids[inactive_links] = - 1
        link_ids[~ inactive_links] = np.arange(
            link_id_offset, link_id_offset + active_link_count,
            dtype=np.int)
    link_ids.shape = (shape[0] - 2, shape[1] - 1)
    return link_ids


def west_links(shape):
    """
    >>> from landlab.utils.structured_grid import west_links
    >>> west_links((3, 4))
    array([[-1,  8,  9, 10],
           [-1, 11, 12, 13],
           [-1, 14, 15, 16]])
    """
    link_ids = horizontal_link_ids(shape)
    link_ids.shape = (shape[0], shape[1] - 1)
    return np.hstack((- np.ones((shape[0], 1), dtype=np.int), link_ids))


def north_links(shape):
    """
    >>> from landlab.utils.structured_grid import north_links
    >>> north_links((3, 4))
    array([[ 0,  1,  2,  3],
           [ 4,  5,  6,  7],
           [-1, -1, -1, -1]])
    """
    link_ids = vertical_link_ids(shape)
    link_ids.shape = (shape[0] - 1, shape[1])
    return np.vstack((link_ids, - np.ones((1, shape[1]), dtype=np.int)))


def south_links(shape):
    """
    >>> from landlab.utils.structured_grid import south_links
    >>> south_links((3, 4))
    array([[-1, -1, -1, -1],
           [ 0,  1,  2,  3],
           [ 4,  5,  6,  7]])
    """
    link_ids = vertical_link_ids(shape)
    link_ids.shape = (shape[0] - 1, shape[1])
    return np.vstack((- np.ones((1, shape[1]), dtype=np.int), link_ids))

def east_links(shape):
    """
    >>> from landlab.utils.structured_grid import east_links
    >>> east_links((3, 4))
    array([[ 8,  9, 10, -1],
           [11, 12, 13, -1],
           [14, 15, 16, -1]])
    """
    link_ids = horizontal_link_ids(shape)
    link_ids.shape = (shape[0], shape[1] - 1)
    return np.hstack((link_ids, - np.ones((shape[0], 1), dtype=int)))


def active_north_links(shape, node_status=None):
    """
    >>> from landlab.utils.structured_grid import active_north_links
    >>> active_north_links((3, 4))
    array([[-1,  0,  1, -1],
           [-1,  2,  3, -1],
           [-1, -1, -1, -1]])
    """
    active_north_links = np.empty(shape, dtype=int)
    try:
        links = vertical_active_link_ids(shape, node_status=node_status)
    except ValueError:
        pass
    links.shape = (shape[0] - 1, shape[1] - 2)
    active_north_links[:-1, 1:-1] = links
    active_north_links[:, (0, -1)] = -1
    active_north_links[-1, :] = -1

    return active_north_links


def active_south_links(shape, node_status=None):
    """
    >>> from landlab.utils.structured_grid import active_south_links
    >>> active_south_links((3, 4))
    array([[-1, -1, -1, -1],
           [-1,  0,  1, -1],
           [-1,  2,  3, -1]])
    """
    active_south_links = np.empty(shape, dtype=int)
    links = vertical_active_link_ids(shape, node_status=node_status)
    links.shape = (shape[0] - 1, shape[1] - 2)
    active_south_links[1:, 1:-1] = links
    active_south_links[:, (0, -1)] = -1
    active_south_links[0, :] = -1

    return active_south_links


def active_west_links(shape, node_status=None):
    """
    >>> from landlab.utils.structured_grid import active_west_links
    >>> active_west_links((3, 4))
    array([[-1, -1, -1, -1],
           [-1,  4,  5,  6],
           [-1, -1, -1, -1]])
    """
    active_west_links = np.empty(shape, dtype=int)
    try:
        active_west_links[1:-1, 1:] = horizontal_active_link_ids(
            shape, node_status=node_status)
    except ValueError:
        pass
    active_west_links[(0, -1), :] = -1
    active_west_links[:, 0] = -1

    return active_west_links


def active_east_links(shape, node_status=None):
    """
    >>> from landlab.utils.structured_grid import active_east_links
    >>> active_east_links((3, 4))
    array([[-1, -1, -1, -1],
           [ 4,  5,  6, -1],
           [-1, -1, -1, -1]])
    """
    active_east_links = np.empty(shape, dtype=int)
    active_east_links.fill(-999)
    try:
        active_east_links[1:-1, :-1] = horizontal_active_link_ids(
            shape, node_status=node_status)
    except ValueError:
        pass
    active_east_links[(0, -1), :] = -1
    active_east_links[:, -1] = -1

    return active_east_links


def outlink_count_per_node(shape):
    link_count = np.empty(shape, dtype=np.int)
    link_count[:-1, :-1] = 2
    link_count[-1, :-1] = 1
    link_count[:-1, -1] = 1
    link_count[-1, -1] = 0
    return np.ravel(link_count)


def inlink_count_per_node(shape):
    link_count = np.empty(shape, dtype=np.int)
    link_count[1:, 1:] = 2
    link_count[0, 1:] = 1
    link_count[1:, 0] = 1
    link_count[0, 0] = 0
    return np.ravel(link_count)


def active_outlink_count_per_node(shape):
    link_count = np.empty(shape, dtype=np.int)
    link_count[1:-1, 1:-1] = 2
    link_count[0, :] = 1
    link_count[-1, :] = 0
    link_count[:, 0] = 1
    link_count[:, -1] = 0

    link_count[0, 0] = 0
    link_count[-1, 0] = 0

    return np.ravel(link_count)


def active_inlink_count_per_node(shape):
    link_count = np.empty(shape, dtype=np.int)
    link_count[1:-1, 1:-1] = 2
    link_count[0, :] = 0
    link_count[-1, :] = 1
    link_count[:, 0] = 0
    link_count[:, -1] = 1

    link_count[0, -1] = 0
    link_count[-1, -1] = 0

    return np.ravel(link_count)


def setup_outlink_matrix(shape, return_count=True):
    links = outlinks(shape)
    if return_count:
        return (links, outlink_count_per_node(shape))
    else:
        return links


def setup_inlink_matrix(shape, return_count=True):
    links = inlinks(shape)
    if return_count:
        return (links, inlink_count_per_node(shape))
    else:
        return links


def setup_active_outlink_matrix(shape, node_status=None, return_count=True):
    links = active_outlinks(shape, node_status=node_status)
    if return_count:
        return links, active_outlink_count_per_node(shape)
    else:
        return links


def setup_active_inlink_matrix(shape, node_status=None, return_count=True):
    """Active links entering nodes.

    Return the IDs of the active links that enter each node of a grid. The shape of the returned array is
    (2, *N*) where *N* is the number of nodes in the grid. The first row contains the link ID entering the
    node from the bottom, and the second row the link entering the node from the left.

    Use the *return_count* keyword to, in addition to the link IDs, return the number of active links attached
    to each grid node.

    Use the *node_status_array* keyword to specify the status for each of the grid's nodes. If not given, each
    of the perimeter nodes is assumed to be `FIXED_VALUE_BOUNDARY`.

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
    Get the active link IDs for a grid of 3 nodes by 4 nodes. The first row list links entering nodes from the
    bottom, and the second links entering from the left.

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


def node_index_with_halo(shape, halo_indices=BAD_INDEX_VALUE):
    """
    >>> from landlab.utils.structured_grid import node_index_with_halo
    >>> node_index_with_halo((2, 3), halo_indices=-1)
    array([[-1, -1, -1, -1, -1],
           [-1,  0,  1,  2, -1],
           [-1,  3,  4,  5, -1],
           [-1, -1, -1, -1, -1]])
    """
    shape_with_halo = np.array(shape) + 2

    ids = np.empty(shape_with_halo, dtype=np.int)

    (interiors, boundaries) = (interior_nodes(shape_with_halo),
                               boundary_nodes(shape_with_halo))

    ids.flat[interiors] = xrange(interior_node_count(shape_with_halo))
    ids.flat[boundaries] = halo_indices

    return ids


def cell_index_with_halo(shape, halo_indices=BAD_INDEX_VALUE,
                        inactive_indices=None):
    """
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
    shape = (ids_with_halo.shape[0] - 2, ids_with_halo.shape[1] - 2)
    kwds = {
        'strides': ids_with_halo.strides,
        'buffer': ids_with_halo,
        'dtype': ids_with_halo.dtype,
    }

    kwds['offset'] = ids_with_halo.itemsize * (ids_with_halo.shape[1])
    west_ids = np.ndarray(shape, **kwds)

    kwds['offset'] = ids_with_halo.itemsize * (ids_with_halo.shape[1] + 2)
    east_ids = np.ndarray(shape, **kwds)

    kwds['offset'] = ids_with_halo.itemsize
    south_ids = np.ndarray(shape, **kwds)

    kwds['offset'] = ids_with_halo.itemsize * (ids_with_halo.shape[1] * 2 + 1)
    north_ids = np.ndarray(shape, **kwds)

    return np.vstack(
        (east_ids.flat, north_ids.flat, west_ids.flat, south_ids.flat))


def _centered_node_ids(ids_with_halo):
    shape = (ids_with_halo.shape[0] - 2, ids_with_halo.shape[1] - 2)
    kwds = {'strides': ids_with_halo.strides,
            'buffer': ids_with_halo,
            'dtype': ids_with_halo.dtype, }

    kwds['offset'] = ids_with_halo.itemsize * (ids_with_halo.shape[1] + 1)
    return np.ndarray(shape, **kwds)


def neighbor_node_ids(shape, inactive=BAD_INDEX_VALUE):
    return linked_neighbor_node_ids(shape, [], inactive=inactive)


def linked_neighbor_node_ids(shape, closed_boundary_nodes,
                             open_boundary_nodes=None,
                             inactive=BAD_INDEX_VALUE):
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
    open_boundary_neighbors = neighbors[:, open_boundary_nodes]
    is_open_boundary_neighbor = _find_open_boundary_neighbors(
        neighbors, open_boundary_nodes)
    nodes = np.choose(is_open_boundary_neighbor, (open_boundary_neighbors,
                                                  value))
    neighbors[:, open_boundary_nodes] = nodes


def _find_open_boundary_neighbors(neighbors, open_boundary_nodes):
    open_boundary_neighbors = neighbors[:, open_boundary_nodes]
    is_open_boundary_neighbor = np.in1d(open_boundary_neighbors,
                                        open_boundary_nodes)
    is_open_boundary_neighbor.shape = (neighbors.shape[0],
                                       len(open_boundary_nodes))
    return is_open_boundary_neighbor


def neighbor_node_array(shape, **kwds):
    """
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
    closed_boundary_nodes = kwds.pop('closed_boundary_nodes', [])
    open_boundary_nodes = kwds.get('open_boundary_nodes', [])

    if len(closed_boundary_nodes) > 0 or len(open_boundary_nodes):
        neighbors = linked_neighbor_node_ids(shape, closed_boundary_nodes,
                                             **kwds)
    else:
        neighbors = neighbor_node_ids(shape, **kwds)

    return neighbors


def neighbor_cell_array(shape, out_of_bounds=BAD_INDEX_VALUE, contiguous=True):
    """
    >>> from landlab.utils.structured_grid import neighbor_cell_array
    >>> neighbors = neighbor_cell_array((2, 3), out_of_bounds=-1)
    >>> neighbors
    array([], dtype=int64)

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

        neighbors = np.vstack((
            ids[1:shape[0] + 1, 2:].flat,
            ids[2:, 1:shape[1] + 1].flat,
            ids[1:shape[0] + 1, :shape[1]].flat,
            ids[:shape[0], 1:shape[1] + 1].flat,)).T
        if contiguous:
            return neighbors.copy()
        else:
            return neighbors
    else:
        return np.array([], dtype=np.int)


def diagonal_node_array(shape, out_of_bounds=BAD_INDEX_VALUE, contiguous=True,
                        boundary_node_mask=None):
    """
    Creates a list of IDs of the diagonal cells to each cell, as a 2D array.
    Only interior cells are assigned neighbors; boundary cells get -1 for
    each neighbor.  The order of the diagonal cells is [topright, topleft,
    bottomleft, bottomright].

    NG didn't touch this, but she thinks this should be nodes, not cells.

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
    ids = node_index_with_halo(shape, halo_indices=out_of_bounds)

    diags = np.vstack((
        ids[2:, 2:].flat,
        ids[2:, :shape[1]].flat,
        ids[:shape[0], :shape[1]].flat,
        ids[:shape[0], 2:].flat,)).T

    if boundary_node_mask is not None:
        boundaries = np.empty(4, dtype=np.int)
        boundaries.fill(boundary_node_mask)
        diags[boundary_nodes(shape)] = boundaries

    if contiguous:
        return diags.copy()
    else:
        return diags


def diagonal_cell_array(shape, out_of_bounds=BAD_INDEX_VALUE, contiguous=True):
    """
    Construct a matrix of cell indices to each diagonally adjacent cell of a
    structured grid. If a cell does not have a diagonal neighbor, set the
    index for that neighbor to *out_of_bounds*.

    An grid without any cells returns an empty array.
    >>> from landlab.utils.structured_grid import diagonal_cell_array
    >>> diags = diagonal_cell_array((2, 3), out_of_bounds=-1)
    >>> diags
    array([], dtype=int64)

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

        diags = np.vstack((
            ids[2:, 2:].flat,
            ids[2:, :shape[1]].flat,
            ids[:shape[0], :shape[1]].flat,
            ids[:shape[0], 2:].flat,)).T
        if contiguous:
            return diags.copy()
        else:
            return diags
    else:
        return np.array([], dtype=np.int)


def diagonal_array_slow(shape, out_of_bounds=BAD_INDEX_VALUE):
    (nrows, ncols) = shape
    ncells = shape[0] * shape[1]
    diagonal_cells = - np.ones([ncells, 4], dtype=np.int)
    for r in xrange( 1, nrows-1 ):
        for c in xrange( 1, ncols-1 ):
            cell_id = r * ncols + c
            diagonal_cells[cell_id,2] = cell_id - ncols - 1 # bottom left
            diagonal_cells[cell_id,0] = cell_id + ncols + 1 # top right
            diagonal_cells[cell_id,3] = cell_id - ncols + 1 # bottom right
            diagonal_cells[cell_id,1] = cell_id + ncols - 1 # top left
    return diagonal_cells


def has_boundary_neighbor(neighbors, diagonals,
                          out_of_bounds=BAD_INDEX_VALUE):
    """
    DEJH thinks this method is broken since terminology update: it returns closed
    neihbors, not boundary neighbors.
    """
    return (out_of_bounds in neighbors |
            out_of_bounds in diagonals)


def has_boundary_neighbor_slow(neighbors, diagonals, out_of_bounds=BAD_INDEX_VALUE):
    # nbr_nodes=self.get_neighbor_list(id)
    # diag_nbrs=self.get_diagonal_list(id)

    i=0
    while i < 4 and neighbors[i] != out_of_bounds:
        i += 1

    if i < 4:
        return True
    else:
        r = 0
        while r < 4 and diagonals[r] != out_of_bounds:
            r += 1

    if r < 4 :
        return True
    else:
        return False


def reshape_array(shape, u, flip_vertically=False, copy=False):
    """
    >>> from landlab.utils.structured_grid import reshape_array
    >>> x = np.arange(12.)
    >>> y = reshape_array((3, 4), x)
    >>> y.shape
    (3, 4)
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
    reshaped_u = u.view()

    try:
        reshaped_u.shape = shape
    except ValueError:
        raise

    if flip_vertically:
        flipped_u = reshaped_u[::-1, :]
        if copy:
            return flipped_u.copy()
        else:
            return flipped_u
    else:
        if copy:
            return reshaped_u.copy()
        else:
            return reshaped_u

def nodes_around_points_on_unit_grid(shape, coords, mode='raise'):
    """
    Returns the nodes around a point on a structured grid with unit spacing
    and zero origin.

    >>> from landlab.utils.structured_grid import nodes_around_points_on_unit_grid
    >>> nodes_around_points_on_unit_grid((3, 3), (.1, .1))
    array([0, 3, 4, 1])

    >>> nodes_around_points_on_unit_grid((3, 3), (1., 1.))
    array([4, 7, 8, 5])
    """
    if isinstance(coords[0], np.ndarray):
        (rows, cols) = (coords[0].astype(np.int), coords[1].astype(np.int))
    else:
        (rows, cols) = (int(coords[0]), int(coords[1]))

    try:
        return np.ravel_multi_index(((rows, rows + 1, rows + 1, rows),
                                     (cols, cols, cols + 1, cols + 1)),
                                    shape, mode=mode).T
    except ValueError:
        raise

def nodes_around_points(shape, coords, spacing=(1., 1.),
                        origin=(0., 0.)):
    """
    Returns the nodes around a point on a structured grid with row and column
    *spacing*, and *origin*.

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
    try:
        return nodes_around_points_on_unit_grid(
            shape,
            ((coords[0] - origin[0]) / spacing[0],
             (coords[1] - origin[1]) / spacing[1]))
    except ValueError:
        raise

def nodes_around_point(shape, coords, spacing=(1., 1.)):
    node_id = int(coords[0] // spacing[0] * shape[1] + coords[1] // spacing[1])
    if node_id + shape[1] + 1 >= shape[0] * shape[1] or node_id < 0:
        raise ValueError('invalid entry in coordinates array')

    return np.array([node_id, node_id + shape[1], node_id + shape[1] + 1,
                     node_id + 1])

def interior_node_id_to_node_id(shape, core_node_ids):
    """
    Converts the id of an interior node ID (i.e., if just the interior nodes
    were numbered) to a node ID.
    """
    IGW = shape[1]-2
    real_ID = (core_node_ids//IGW + 1) * shape[1] + (core_node_ids%IGW) + 1
    assert np.all(real_ID < shape[0]*shape[1])
    return real_ID.astype(int)

def node_id_to_interior_node_id(shape, node_ids):
    """
    Converts a node ID to the id of an interior node ID (i.e., if just the
    interior nodes were numbered)
    """
    ncols = shape[1]
    interior_ID = (node_ids//ncols - 1)*(ncols-2) + (node_ids%ncols) - 1
    if np.any(interior_ID < 0) or np.any(interior_ID >= (shape[0]-2)*(shape[1]-2)):
        print "One of the supplied nodes was outside the interior grid!"
        raise NameError()
    else:
        return interior_ID.astype(int)

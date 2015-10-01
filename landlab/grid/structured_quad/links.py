import numpy as np
from . import nodes
from ..base import CORE_NODE, FIXED_GRADIENT_BOUNDARY, FIXED_VALUE_BOUNDARY
from ..unstructured.links import LinkGrid
from ...core.utils import as_id_array


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
    array([[0, 1, 2, 3],
           [4, 5, 6, 7]])
    """
    link_ids = np.arange(number_of_vertical_links(shape), dtype=np.int)
    return link_ids.reshape(shape_of_vertical_links(shape))


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
    array([[ 8,  9, 10],
           [11, 12, 13],
           [14, 15, 16]])
    """
    link_ids = (np.arange(number_of_horizontal_links(shape), dtype=np.int) +
                number_of_vertical_links(shape))
    return link_ids.reshape(shape_of_horizontal_links(shape))


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
    array([[ 0,  1,  2,  3],
           [ 4,  5,  6,  7],
           [-1, -1, -1, -1]])
    >>> horiz
    array([[ 8,  9, 10, -1],
           [11, 12, 13, -1],
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
           [ 0,  1,  2,  3],
           [ 4,  5,  6,  7]])
    >>> horiz
    array([[-1,  8,  9, 10],
           [-1, 11, 12, 13],
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
    array([ 8,  9, 10,  0,  1, 11,  2, 12,  3, 13,  4,  5, 14,  6, 15,  7, 16])
    >>> offset
    array([ 0,  0,  1,  2,  3,  4,  6,  8, 10, 11, 13, 15, 17])

    The links entering the 1st, 5th, and last node. The first node does not
    have any links entering it.

    >>> offset[0] == offset[1]
    True
    >>> for link in [4, 11]: links[offset[link]:offset[link + 1]]
    array([0])
    array([ 7, 16])
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
    array([ 0,  8,  1,  9,  2, 10,  3,  4, 11,  5, 12,  6, 13,  7, 14, 15, 16])
    >>> offset
    array([ 0,  2,  4,  6,  7,  9, 11, 13, 14, 15, 16, 17, 17])

    The links leaving the 1st, 8th, and last node. The last node does not have
    any links leaving it.

    >>> offset[11] == offset[12]
    True
    >>> for link in [0, 7]: links[offset[link]:offset[link + 1]]
    array([0, 8])
    array([7])
    """
    (out_vert, out_horiz) = _node_out_link_ids(shape)
    _node_link_ids = np.vstack((out_vert.flat, out_horiz.flat)).T
    offset = np.empty(nodes.number_of_nodes(shape) + 1, dtype=int)
    np.cumsum(number_of_out_links_per_node(shape), out=offset[1:])
    offset[0] = 0
    return _node_link_ids[_node_link_ids >= 0], offset


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
    array([ 0,  8,  8,  1,  9,  9,  2, 10, 10,  3,  0,  4, 11,  1, 11,  5, 12,
            2, 12,  6, 13,  3, 13,  7,  4, 14,  5, 14, 15,  6, 15, 16,  7, 16])
    >>> offset
    array([ 0,  2,  5,  8, 10, 13, 17, 21, 24, 26, 29, 32, 34])

    The links attached to node 0

    >>> links[offset[0]:offset[1]]
    array([0, 8])

    The links attached to node 5

    >>> links[offset[5]:offset[6]]
    array([ 1, 11,  5, 12])
    """
    (in_vert, in_horiz) = _node_in_link_ids(shape)
    (out_vert, out_horiz) = _node_out_link_ids(shape)
    _node_link_ids = np.vstack((in_vert.flat, in_horiz.flat,
                                out_vert.flat, out_horiz.flat)).T

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
    >>> node_id_at_link_start((3, 4))
    array([ 0,  1,  2,  3,  4,  5,  6,  7,  0,  1,  2,  4,  5,  6,  8,  9, 10])
    """
    all_node_ids = nodes.node_ids(shape)
    return np.concatenate((all_node_ids[:-1, :].flat,
                           all_node_ids[:, :-1].flat))


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
    >>> node_id_at_link_end((3, 4))
    array([ 4,  5,  6,  7,  8,  9, 10, 11,  1,  2,  3,  5,  6,  7,  9, 10, 11])
    """
    all_node_ids = nodes.node_ids(shape)
    return np.concatenate((all_node_ids[1:, :].flat, all_node_ids[:, 1:].flat))


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
    >>> is_active_link((3, 4), status)
    array([False, False, False, False, False, False, False, False, False,
           False, False, False,  True, False, False, False, False], dtype=bool)
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
    """returns link ids of active links.

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
    >>> rmg = RasterModelGrid(3, 4)
    >>> rmg.set_closed_boundaries_at_grid_edges(True, True, True, True)
    >>> status = rmg.status_at_node
    >>> active_link_ids((3, 4), status)
    array([12])
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
    >>> rmg = RasterModelGrid(4, 5)
    >>> z = np.arange(0, rmg.number_of_nodes)
    >>> s = np.arange(0, rmg.number_of_links)
    >>> rmg['node']['topographic__elevation'] = z
    >>> rmg['link']['topographic__slope'] = s
    >>> rmg.set_fixed_link_boundaries_at_grid_edges(True, True, True, True)
    >>> is_fixed_link(rmg.shape, rmg.status_at_node)
    array([False,  True,  True,  True, False, False, False, False, False,
           False, False,  True,  True,  True, False, False, False, False,
           False,  True, False, False,  True,  True, False, False,  True,
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
    >>> rmg['node']['topographic__elevation'] = z
    >>> rmg['link']['topographic__slope'] = s
    >>> rmg.set_fixed_link_boundaries_at_grid_edges(True, True, True, True)
    >>> fixed_link_ids(rmg.shape, rmg.status_at_node)
    array([ 1,  2,  3, 11, 12, 13, 19, 22, 23, 26])
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
    >>> active_ids = active_link_ids((4,5), status)
    >>> horizontal_active_link_ids((4,5), active_ids)
    array([-1, -1, -1, -1, -1, 20, 21, -1, -1, 24, 25, -1, -1, -1, -1, -1])
    """
    # For horizontal links, we need to start with a list of '-1' indices with
    # length of number_of_links
    horizontal_links = np.ones(number_of_links(shape)) * bad_index_value

    # We will need the number of rows and columns from input argument 'shape'
    rows = shape[0]
    cols = shape[1]

    # In a structured quad, the minimum horizontal link id is equal to the
    # number of columns * (number of rows - 1)
    min_hori_id = cols * (rows - 1)

    # We will use list comprehension to get *just* the horizontal link ids
    # from the active_link_ids input argument.
    horizontal_ids = [i for i in active_ids if i >= min_hori_id]

    # In the array of '-1' we input the horizontal active link ids
    horizontal_links[horizontal_ids] = horizontal_ids

    # To get an array of len number_of_horizontal_links, we need to clip off
    # the number of vertical links. We do this by starting at the "minimum
    # horizontal link id" found above and going to the end of the list.
    horizontal_links = horizontal_links[min_hori_id:]
    horizontal_links = horizontal_links.astype(int)

    # Return an array with length of number_of_vertical_links that has '-1' for
    # inactive/fixed links and the active link id for active links
    return horizontal_links


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
          *--23-->o------>o------>o--26-->*
          ^       ^       ^       ^       ^
          I       V       V       V       I
          |       |       |       |       |
          *--19-->o------>o------>o--22-->*
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
    >>> rmg['node']['topographic__elevation'] = numpy.arange(
    ...     0, rmg.number_of_nodes)
    >>> rmg['link']['topographic__slope'] = numpy.arange(
    ...     0, rmg.number_of_links)
    >>> rmg.set_fixed_link_boundaries_at_grid_edges(True, True, True, True)
    >>> status = rmg.status_at_node
    >>> fixed_ids = fixed_link_ids((4,5), status)
    >>> horizontal_fixed_link_ids((4,5), fixed_ids)
    array([-1, -1, -1, -1, 19, -1, -1, 22, 23, -1, -1, 26, -1, -1, -1, -1])
    """
    # For horizontal links, we need to start with a list of '-1' indices with
    # length of number_of_links
    horizontal_links = np.ones(number_of_links(shape)) * bad_index_value

    # We will need the number of rows and columns from input argument 'shape'
    rows = shape[0]
    cols = shape[1]

    # In a structured quad, the minimum horizontal link id is equal to the
    # number of columns * (number of rows - 1)
    min_hori_id = cols * (rows - 1)

    # We will use list comprehension to get *just* the horizontal link ids
    # from the fixed_link_ids input argument.
    horizontal_ids = [i for i in fixed_ids if i >= min_hori_id]

    # In the array of '-1' we input the horizontal fixed link ids
    horizontal_links[horizontal_ids] = horizontal_ids

    # To get an array of len number_of_horizontal_links, we need to clip off
    # the number of vertical links. We do this by starting at the "minimum
    # horizontal link id" found above and going to the end of the list.
    horizontal_links = horizontal_links[min_hori_id:]

    # Return an array with length of number_of_vertical_links that has '-1' for
    # inactive links and the active link id for active links
    return as_id_array(horizontal_links)


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
    >>> rmg = RasterModelGrid(4, 5)
    >>> rmg.set_closed_boundaries_at_grid_edges(True, True, True, True)
    >>> status = rmg.status_at_node
    >>> active_ids = active_link_ids((4,5), status)
    >>> vertical_active_link_ids((4,5), active_ids)
    array([-1, -1, -1, -1, -1, -1,  6,  7,  8, -1, -1, -1, -1, -1, -1])
    """
    # Set up an array of '-1' indices with length of number_of_vertical_links
    vertical_links = np.ones(number_of_vertical_links(shape)) * bad_index_value

    # We will need the number of rows and columns from input argument 'shape'
    rows = shape[0]
    cols = shape[1]

    # In a structured quad, the maximum vertical link id is one less than the
    # number of columns * (number of rows - 1)
    max_vert_id = cols * (rows - 1)

    # We will use list comprehension to get *just* the vertical link ids
    # from the active_link_ids input argument.
    vertical_ids = [i for i in active_ids if i < max_vert_id]

    # In the array of '-1's, we input the active link ids.
    vertical_links[vertical_ids] = vertical_ids
    vertical_links = vertical_links.astype(int)

    # Return an array with length of number_of_vertical_links that has '-1' for
    # inactive links and the active link id for active links
    return vertical_links


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
        I      11      12      13       I
        |       |       |       |       |
        *---H-->o---H-->o---H-->o---H-->*
        ^       ^       ^       ^       ^
        I       V       V       V       I
        |       |       |       |       |
        *---H-->o---H-->o---H-->o---H-->*
        ^       ^       ^       ^       ^
        I       1       2       3       I
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
    >>> rmg = RasterModelGrid(4, 5)
    >>> rmg['node']['topographic__elevation'] = numpy.arange(
    ...     0, rmg.number_of_nodes)
    >>> rmg['link']['topographic__slope'] = numpy.arange(
    ...     0, rmg.number_of_links)
    >>> rmg.set_fixed_link_boundaries_at_grid_edges(True, True, True, True)
    >>> status = rmg.status_at_node
    >>>
    >>> fixed_ids = fixed_link_ids((4,5), status)
    >>> vertical_fixed_link_ids((4,5), fixed_ids)
    array([-1,  1,  2,  3, -1, -1, -1, -1, -1, -1, -1, 11, 12, 13, -1])
    """
    # Set up an array of '-1' indices with length of number_of_vertical_links
    vertical_links = np.ones(number_of_vertical_links(shape)) * bad_index_value

    # We will need the number of rows and columns from input argument 'shape'
    rows = shape[0]
    cols = shape[1]

    # In a structured quad, the maximum vertical link id is one less than the
    # number of columns * (number of rows - 1)
    max_vert_id = cols * (rows - 1)

    # We will use list comprehension to get *just* the vertical link ids
    # from the active_link_ids input argument.
    vertical_ids = [i for i in fixed_ids if i < max_vert_id]

    # In the array of '-1's, we input the active link ids.
    vertical_links[vertical_ids] = vertical_ids

    # Return an array with length of number_of_vertical_links that has '-1' for
    # inactive links and the active link id for active links

    return as_id_array(vertical_links)


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



        *--23-->*--24-->*--25-->*--26-->*



        *--19-->*--20-->*--21-->*--22-->*



        *--15-->*--16-->*--17-->*--18-->*

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
    array([-1, -1, -1, -1, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26])
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



        *--23-->*--24-->*--25-->*--26-->*



        *--19-->*--20-->*--21-->*--22-->*



        *--15-->*--16-->*--17-->*--18-->*

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
    array([-1, 15, 16, 17, -1, 19, 20, 21, -1, 23, 24, 25, -1, 27, 28, 29])
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



        *--23-->*--24-->*--25-->*--26-->*



        *--19-->*--20-->*--21-->*--22-->*



        *--15-->*--16-->*--17-->*--18-->*

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
    array([19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, -1, -1, -1, -1])
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



        *--23-->*--24-->*--25-->*--26-->*



        *--19-->*--20-->*--21-->*--22-->*



        *--15-->*--16-->*--17-->*--18-->*

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
    array([16, 17, 18, -1, 20, 21, 22, -1, 24, 25, 26, -1, 28, 29, 30, -1])
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
        [S, W, N, E].

    Examples
    --------
    Sample grid, giving neighbors for link ID 20::

        *------>*------>*------>*------>*



        *------>*--24-->*------>*------>*



        *--19-->*--20-->*--21-->*------>*



        *------>*--16-->*------>*------>*

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
    array([[-1, -1, 19, 16],
           [-1, 15, 20, 17],
           [-1, 16, 21, 18],
           [-1, 17, 22, -1],
           [15, -1, 23, 20],
           [16, 19, 24, 21],
           [17, 20, 25, 22],
           [18, 21, 26, -1],
           [19, -1, 27, 24],
           [20, 23, 28, 25],
           [21, 24, 29, 26],
           [22, 25, 30, -1],
           [23, -1, -1, 28],
           [24, 27, -1, 29],
           [25, 28, -1, 30],
           [26, 29, -1, -1]])
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
    neighbor_array = np.array([south, west, north, east])

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
    ndarray :
        Array of 4 horizontal link neighbors for a given link ID. Returned in
        [S, W, N, E]. Returns array for only ACTIVE horizontal links.

    Examples
    --------
    Sample grid, giving neighbors for link ID 20::

        *------>*------>*------>*------>*



        *------>*--24-->*------>*------>*



        *------>*--20-->*--21-->*------>*



        *------>*------>*------>*------>*


    .. note::

        Only horizontal links are shown. When no neighbor is found,
        bad_index_value is returned.

        ``*`` indicates nodes

        Numeric values correspond to the horizontal :any:`ACTIVE_LINK` IDs.

    >>> from landlab import RasterModelGrid
    >>> from landlab.grid.structured_quad.links import *
    >>> rmg = RasterModelGrid(4, 5)
    >>> rmg.set_closed_boundaries_at_grid_edges(True, True, True, True)
    >>> active_ids = active_link_ids(rmg.shape, rmg.status_at_node)
    >>> horizontal_ids = horizontal_active_link_ids(
    ...     rmg.shape, active_ids)
    >>> d4_horizontal_active_link_neighbors(rmg.shape, horizontal_ids)
    array([[-1, -1, 24, 21],
           [-1, 20, 25, -1],
           [20, -1, -1, 25],
           [21, 24, -1, -1]])
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
       10       11      12      13      14
        |       |       |       |       |
        *       *       *       *       *
        ^       ^       ^       ^       ^
        5       6       7       8       9
        |       |       |       |       |
        *       *       *       *       *
        ^       ^       ^       ^       ^
        0       1       2       3       4
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
    >>> vertical_south_link_neighbor(rmg.shape, vertical_links)
    array([-1, -1, -1, -1, -1,  0,  1,  2,  3,  4,  5,  6,  7,  8,  9])
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
       10       11      12      13      14
        |       |       |       |       |
        *       *       *       *       *
        ^       ^       ^       ^       ^
        5       6       7       8       9
        |       |       |       |       |
        *       *       *       *       *
        ^       ^       ^       ^       ^
        0       1       2       3       4
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
    array([-1,  0,  1,  2,  3, -1,  5,  6,  7,  8, -1, 10, 11, 12, 13])
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
       10       11      12      13      14
        |       |       |       |       |
        *       *       *       *       *
        ^       ^       ^       ^       ^
        5       6       7       8       9
        |       |       |       |       |
        *       *       *       *       *
        ^       ^       ^       ^       ^
        0       1       2       3       4
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
    array([ 5,  6,  7,  8,  9, 10, 11, 12, 13, 14, -1, -1, -1, -1, -1])

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
       10       11      12      13      14
        |       |       |       |       |
        *       *       *       *       *
        ^       ^       ^       ^       ^
        5       6       7       8       9
        |       |       |       |       |
        *       *       *       *       *
        ^       ^       ^       ^       ^
        0       1       2       3       4
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
    array([ 1,  2,  3,  4, -1,  6,  7,  8,  9, -1, 11, 12, 13, 14, -1])
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
        [S, W, N, E].

    Examples
    --------
    The following example uses this grid::

        *       *       *       *       *
        ^       ^       ^       ^       ^
       10       11      12      13      14
        |       |       |       |       |
        *       *       *       *       *
        ^       ^       ^       ^       ^
        5       6       7       8       9
        |       |       |       |       |
        *       *       *       *       *
        ^       ^       ^       ^       ^
        0       1       2       3       4
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
    array([[-1, -1,  5,  1],
           [-1,  0,  6,  2],
           [-1,  1,  7,  3],
           [-1,  2,  8,  4],
           [-1,  3,  9, -1],
           [ 0, -1, 10,  6],
           [ 1,  5, 11,  7],
           [ 2,  6, 12,  8],
           [ 3,  7, 13,  9],
           [ 4,  8, 14, -1],
           [ 5, -1, -1, 11],
           [ 6, 10, -1, 12],
           [ 7, 11, -1, 13],
           [ 8, 12, -1, 14],
           [ 9, 13, -1, -1]])
    """
    south = vertical_south_link_neighbor(shape, vertical_ids, bad_index_value)
    west = vertical_west_link_neighbor(shape, vertical_ids, bad_index_value)
    north = vertical_north_link_neighbor(shape, vertical_ids, bad_index_value)
    east = vertical_east_link_neighbor(shape, vertical_ids, bad_index_value)
    neighbor_array = np.array([south, west, north, east])
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
        Returned in [S, W, N, E].

    Examples
    --------
    The following example uses this grid::

        *       *       *       *       *
        ^       ^       ^       ^       ^
       10       11      12      13      14
        |       |       |       |       |
        *       *       *       *       *
        ^       ^       ^       ^       ^
        5       6       7       8       9
        |       |       |       |       |
        *       *       *       *       *
        ^       ^       ^       ^       ^
        0       1       2       3       4
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
    array([[-1, -1,  6,  2],
           [-1,  1,  7,  3],
           [-1,  2,  8, -1],
           [ 1, -1, 11,  7],
           [ 2,  6, 12,  8],
           [ 3,  7, 13, -1],
           [ 6, -1, -1, 12],
           [ 7, 11, -1, 13],
           [ 8, 12, -1, -1]])
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

        *--28-->*--29-->*--30-->*--31-->*



        *--24-->*--25-->*--26-->*--27-->*



        *--20-->*--21-->*--22-->*--23-->*



        *--15-->*--16-->*--17-->*--18-->*

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
    array([15, 16, 17, 18])
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

        *--28-->*--29-->*--30-->*--31-->*



        *--24-->*--25-->*--26-->*--27-->*



        *--20-->*--21-->*--22-->*--23-->*



        *--15-->*--16-->*--17-->*--18-->*

    .. note::

        Only horizontal links are shown.

        ``*`` indicates nodes

        Numeric values correspond to the horizontal IDs.

    >>> from landlab import RasterModelGrid
    >>> from landlab.grid.structured_quad.links import left_edge_horizontal_ids
    >>> rmg = RasterModelGrid(4, 5)
    >>> shape = rmg.shape
    >>> left_edge_horizontal_ids(shape)
    array([15, 19, 23, 27])
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

        *--28-->*--29-->*--30-->*--31-->*



        *--24-->*--25-->*--26-->*--27-->*



        *--20-->*--21-->*--22-->*--23-->*



        *--15-->*--16-->*--17-->*--18-->*

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

        *--28-->*--29-->*--30-->*--31-->*



        *--24-->*--25-->*--26-->*--27-->*



        *--20-->*--21-->*--22-->*--23-->*



        *--15-->*--16-->*--17-->*--18-->*

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
    array([18, 22, 26, 30])
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
       10       11      12      13      14
        |       |       |       |       |
        *       *       *       *       *
        ^       ^       ^       ^       ^
        5       6       7       8       9
        |       |       |       |       |
        *       *       *       *       *
        ^       ^       ^       ^       ^
        0       1       2       3       4
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
    array([0, 1, 2, 3, 4])
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
       10       11      12      13      14
        |       |       |       |       |
        *       *       *       *       *
        ^       ^       ^       ^       ^
        5       6       7       8       9
        |       |       |       |       |
        *       *       *       *       *
        ^       ^       ^       ^       ^
        0       1       2       3       4
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
    array([ 0,  5, 10])
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
       10       11      12      13      14
        |       |       |       |       |
        *       *       *       *       *
        ^       ^       ^       ^       ^
        5       6       7       8       9
        |       |       |       |       |
        *       *       *       *       *
        ^       ^       ^       ^       ^
        0       1       2       3       4
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
    array([10, 11, 12, 13, 14])
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
       10       11      12      13      14
        |       |       |       |       |
        *       *       *       *       *
        ^       ^       ^       ^       ^
        5       6       7       8       9
        |       |       |       |       |
        *       *       *       *       *
        ^       ^       ^       ^       ^
        0       1       2       3       4
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
    array([ 4,  9, 14])
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

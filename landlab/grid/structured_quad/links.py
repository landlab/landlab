import numpy as np


from . import nodes
from ..base import CORE_NODE, CLOSED_BOUNDARY
from ..unstructured.links import LinkGrid


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
    return (shape[0] - 1,  shape[1])


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
    >>> from landlab.grid.structured_quad.links import shape_of_horizontal_links
    >>> shape_of_horizontal_links((3, 4))
    (3, 3)
    """
    return (shape[0],  shape[1] - 1)


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
    >>> from landlab.grid.structured_quad.links import number_of_horizontal_links
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
    """IDs of vertical links in a structured quad grid.

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
    link_ids = np.arange(number_of_vertical_links(shape))
    return link_ids.reshape(shape_of_vertical_links(shape))


def horizontal_link_ids(shape):
    """IDs of horizontal links in a structured quad grid.

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
    link_ids = (np.arange(number_of_horizontal_links(shape)) +
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
    >>> from landlab.grid.structured_quad.links import (number_of_links_per_node,
    ...                                            number_of_in_links_per_node,
    ...                                            number_of_out_links_per_node)
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
    >>> from landlab.grid.structured_quad.links import number_of_in_links_per_node
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
    >>> from landlab.grid.structured_quad.links import number_of_out_links_per_node
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
    """Link IDs for links leaving each node.

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
    """Link IDs for links entering each node.

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
    """Link IDs for links entering each node.

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

    The links entering the 1st, 5th, and last node,

    >>> for link in [0, 4, 11]: links[offset[link]:offset[link + 1]]
    array([], dtype=int64)
    array([0])
    array([ 7, 16])
    """
    (in_vert, in_horiz) = _node_in_link_ids(shape)
    node_link_ids = np.vstack((in_vert.flat, in_horiz.flat)).T
    #offset = np.cumsum(number_of_in_links_per_node(shape))

    offset = np.empty(nodes.number_of_nodes(shape) + 1, dtype=int)
    np.cumsum(number_of_in_links_per_node(shape), out=offset[1:])
    offset[0] = 0
    return node_link_ids[node_link_ids >= 0], offset


def node_out_link_ids(shape):
    """Link IDs for links leaving each node.

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

    The links leaving the 1st, 8th, and last node,

    >>> for link in [0, 7, 11]: links[offset[link]:offset[link + 1]]
    array([0, 8])
    array([7])
    array([], dtype=int64)
    """
    (out_vert, out_horiz) = _node_out_link_ids(shape)
    node_link_ids = np.vstack((out_vert.flat, out_horiz.flat)).T
    offset = np.empty(nodes.number_of_nodes(shape) + 1, dtype=int)
    np.cumsum(number_of_out_links_per_node(shape), out=offset[1:])
    offset[0] = 0
    return node_link_ids[node_link_ids >= 0], offset


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
    node_link_ids = np.vstack((in_vert.flat, in_horiz.flat, out_vert.flat, out_horiz.flat)).T

    offset = np.empty(nodes.number_of_nodes(shape) + 1, dtype=int)
    np.cumsum(number_of_links_per_node(shape), out=offset[1:])
    offset[0] = 0

    return node_link_ids[node_link_ids >= 0], offset


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
    return np.concatenate((all_node_ids[:-1, :].flat, all_node_ids[:, :-1].flat))


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
    """IDs of active links.

    Parameters
    ----------
    shape : tuple of int
        Shape of grid of nodes.

    Returns
    -------
    ndarray :
        Links IDs at the active links.

    Examples
    --------
    >>> from landlab.grid.structured_quad.nodes import status_with_perimeter_as_boundary
    >>> from landlab.grid.structured_quad.links import is_active_link
    >>> status = status_with_perimeter_as_boundary((3, 4))
    >>> is_active_link((3, 4), status)
    array([False,  True,  True, False, False,  True,  True, False, False,
           False, False,  True,  True,  True, False, False, False], dtype=bool)
    """
    status_at_link_start = node_status.flat[node_id_at_link_start(shape)]
    status_at_link_end = node_status.flat[node_id_at_link_end(shape)]

    return (((status_at_link_start == CORE_NODE) &
             ~ (status_at_link_end == CLOSED_BOUNDARY)) |
            ((status_at_link_end == CORE_NODE) &
             ~ (status_at_link_start == CLOSED_BOUNDARY)))


def active_link_ids(shape, node_status):
    """IDs of active links.

    Parameters
    ----------
    shape : tuple of int
        Shape of grid of nodes.

    Returns
    -------
    ndarray :
        Links IDs at the active links.

    Examples
    --------
    >>> from landlab.grid.structured_quad.nodes import status_with_perimeter_as_boundary
    >>> from landlab.grid.structured_quad.links import active_link_ids
    >>> status = status_with_perimeter_as_boundary((3, 4))
    >>> active_link_ids((3, 4), status)
    array([ 1,  2,  5,  6, 11, 12, 13])
    """
    return np.where(is_active_link(shape, node_status))[0]


class StructuredQuadLinkGrid(LinkGrid):
    def __init__(self, shape):
        link_ends = (node_id_at_link_start(shape), node_id_at_link_end(shape))
        number_of_nodes = np.prod(shape)
        LinkGrid.__init__(self, link_ends, number_of_nodes)

import numpy as np


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
    >>> (vert, horiz) = node_out_link_ids((3, 4))
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
    >>> (vert, horiz) = node_in_link_ids((3, 4))
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
    >>> (links, offset) = node_link_ids((3, 4))
    >>> links
    array([ 0,  8,  8,  1,  9,  9,  2, 10, 10,  3,  0,  4, 11,  1, 11,  5, 12,
            2, 12,  6, 13,  3, 13,  7,  4, 14,  5, 14, 15,  6, 15, 16,  7, 16])
    >>> offset
    array([ 2,  5,  8, 10, 13, 17, 21, 24, 26, 29, 32, 34])

    The links attached to node 0

    >>> links[:offset[0]]
    array([0, 8])

    The links attached to node 5

    >>> links[offset[4]:offset[5]]
    array([ 1, 11,  5, 12])
    """
    (in_vert, in_horiz) = node_in_link_ids(shape)
    (out_vert, out_horiz) = node_out_link_ids(shape)
    node_link_ids = np.vstack((in_vert.flat, in_horiz.flat, out_vert.flat, out_horiz.flat)).T
    offset = np.cumsum(number_of_links_per_node(shape))
    return node_link_ids[node_link_ids >= 0], offset

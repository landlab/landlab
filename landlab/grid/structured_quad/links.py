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
    link_ids = np.arange(number_of_vertical_links(shape), dtype=np.int)
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

    The links entering the 1st, 5th, and last node. The first node does not
    have any links entering it.

    >>> offset[0] == offset[1]
    True
    >>> for link in [4, 11]: links[offset[link]:offset[link + 1]]
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

    The links leaving the 1st, 8th, and last node. The last node does not have
    any links leaving it.

    >>> offset[11] == offset[12]
    True
    >>> for link in [0, 7]: links[offset[link]:offset[link + 1]]
    array([0, 8])
    array([7])
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
    node_status : array_link
        Status of nodes in grid.

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
    if np.prod(shape) != node_status.size:
        raise ValueError('node status array does not match size of grid '
                         '(%d != %d)' % (np.prod(shape), len(node_status)))

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
    node_status : array_link
        Status of nodes in grid.

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
    return np.where(is_active_link(shape, node_status))[0].astype(np.int,
                                                                  copy=False)

    
def horizontal_active_link_ids(shape, active_link_ids, BAD_INDEX_VALUE=-1):
    """Get IDs of horizontal active links.

    Parameters
    ----------
    shape : tuple of int
        Shape of grid of nodes.
    active_link_ids : array of int
        Array of all active link ids
    BAD_INDEX_VALUE: int, optional
        Value assigned to inactive indicies in the array.

    Returns
    -------
    ndarray :
        Link IDs at the HORIZONTAL active links. Length of number_of_horizontal_links.

    Examples
    --------
    >>> from landlab import RasterModelGrid
    >>> from landlab.grid.structured_quad.links import active_link_ids, horizontal_active_link_ids
    >>>
    >>> rmg = RasterModelGrid(4, 5)
    >>> rmg.set_closed_boundaries_at_grid_edges(True, True, True, True)
    >>> status = rmg.node_status 
    >>> 
    >>> active_link_ids = active_link_ids((4,5), status)
    >>> horizontal_active_link_ids((4,5), active_link_ids)
    array([-1, -1, -1, -1, -1, 20, 21, -1, -1, 24, 25, -1, -1, -1, -1, -1])

    

    Example grid: Indicies are given for active horizontal links in the 4x5 grid space.
                "X" indicates inactive links, "V" indicates vertical ids

    ::

          *---X-->*---X-->*---X-->*---X-->*
          ^       ^       ^       ^       ^   
          X       X       X       X       X      
          |       |       |       |       |         
          *---X-->*--24-->*--25-->*---X-->*
          ^       ^       ^       ^       ^       
          X       V       V       V       X       
          |       |       |       |       |       
          *---X-->*--20-->*--21-->*---X-->*
          ^       ^       ^       ^       ^   
          X       X       X       X       X       
          |       |       |       |       |       
          *---X-->*---X-->*---X-->*---X-->*
    """
    # For horizontal links, we need to start with a list of '-1' indices with 
    # length of number_of_links
    horizontal_links = np.ones(number_of_links(shape))*BAD_INDEX_VALUE
    
    # We will need the number of rows and columns from input argument 'shape'
    rows = shape[0]
    cols = shape[1]
    
    # In a structured quad, the minimum horizontal link id is equal to the 
    # number of columns * (number of rows - 1)
    min_hori_id = cols*(rows-1)
    
    # We will use list comprehension to get *just* the horizontal link ids
    # from the active_link_ids input argument.
    horizontal_ids = [i for i in active_link_ids if i>= min_hori_id]
    
    # In the array of '-1' we input the horizontal active link ids
    horizontal_links[horizontal_ids] = horizontal_ids
    
    # To get an array of len number_of_horizontal_links, we need to clip off the 
    # number of vertical links. We do this by starting at the "minimum horizontal
    # link id" found above and going to the end of the list. 
    horizontal_links = horizontal_links[min_hori_id:]
    horizontal_links = horizontal_links.astype(int)
    
    # Return an array with length of number_of_vertical_links that has '-1' for
    # inactive links and the active link id for active links
    return horizontal_links


def vertical_active_link_ids(shape, active_link_ids, BAD_INDEX_VALUE=-1):
    """Get IDs of vertical active links.

    Parameters
    ----------
    shape : tuple of int
        Shape of grid of nodes.
    active_link_ids : array of int
        Array of all active link ids
    BAD_INDEX_VALUE: int, optional
        Value assigned to inactive indicies in the array.

    Returns
    -------
    ndarray :
        Link IDs at the VERTICAL active links. Length of number_of_vertical_links.

    Examples
    --------
    >>> from landlab import RasterModelGrid
    >>> from landlab.grid.structured_quad.links import active_link_ids, vertical_active_link_ids
    >>>
    >>> rmg = RasterModelGrid(4, 5)
    >>> rmg.set_closed_boundaries_at_grid_edges(True, True, True, True)
    >>> status = rmg.node_status 
    >>> 
    >>> active_link_ids = active_link_ids((4,5), status)
    >>> vertical_active_link_ids((4,5), active_link_ids)
    array([-1, -1, -1, -1, -1, -1,  6,  7,  8, -1, -1, -1, -1, -1, -1])
    
    
    Example grid: Indicies are given for active vertical links in the 4x5 grid space.
                "X" indicates inactive links, "H" indicates horizontal ids
                
    ::

          *---X-->*---X-->*---X-->*---X-->*
          ^       ^       ^       ^       ^   
          X       X       X       X       X      
          |       |       |       |       |         
          *---X-->*---H-->*---H-->*---X-->*
          ^       ^       ^       ^       ^       
          X       6       7       8       X       
          |       |       |       |       |       
          *---X-->*---H-->*---H-->*---X-->*
          ^       ^       ^       ^       ^   
          X       X       X       X       X       
          |       |       |       |       |       
          *---X-->*---X-->*---X-->*---X-->*
    """
    # Set up an array of '-1' indices with length of number_of_vertical_links
    vertical_links = np.ones(number_of_vertical_links(shape))*BAD_INDEX_VALUE
    
    # We will need the number of rows and columns from input argument 'shape'
    rows = shape[0]
    cols = shape[1]
    
    # In a structured quad, the maximum vertical link id is one less than the 
    # number of columns * (number of rows - 1)
    max_vert_id = cols*(rows-1)
    
    # We will use list comprehension to get *just* the vertical link ids
    # from the active_link_ids input argument.
    vertical_ids = [i for i in active_link_ids if i < max_vert_id]
    
    # In the array of '-1's, we input the active link ids. 
    vertical_links[vertical_ids] = vertical_ids
    vertical_links = vertical_links.astype(int)

    # Return an array with length of number_of_vertical_links that has '-1' for
    # inactive links and the active link id for active links
    return vertical_links
    
def find_horizontal_south_neighbor(shape, horizontal_link_ids, BAD_INDEX_VALUE=-1):
    """Get IDs of SOUTH, horizontal link neighbor

    Parameters
    ----------
    shape : tuple of int
        Shape of grid of nodes.
    horizontal_active_ids : array of int
        Array of all horizontal link ids - MUST BE ARRAY OF LEN(HORIZONTAL_LINKS)
    BAD_INDEX_VALUE: int, optional
        Value assigned to inactive indicies in the array.

    Returns
    -------
    ndarray :
        Link IDs of south horizontal neighbor active links. Length of number_of_horizontal_links.

    Examples
    --------
    >>> from landlab import RasterModelGrid
    >>> from landlab.grid.structured_quad.links import active_link_ids, horizontal_active_link_ids, find_horizontal_south_neighbor

    >>> rmg = RasterModelGrid(4, 5)

    >>> active_link_ids = active_link_ids((4,5), rmg.node_status)
    >>> horizontal_active_ids = horizontal_active_link_ids((4,5), active_link_ids)
    >>> find_horizontal_south_neighbor((4,5), horizontal_active_ids)
    array([-1, -1, -1, -1, -1, -1, -1, -1, 19, 20, 21, 22, 23, 24, 25, 26])

    
    Example grid: Indicies are given for active horizontal links in the 4x5 grid space.
                  Only horizontal links are shown, asterisks (*) represent nodes.
                  In this example, link 24 has one south neighbor, link 20. In the 
                  find_horizontal_south_neighbor array, the index of 20 (the south link)
                  corresponds with link 24. Similarly, for link 25, the south neighbor
                  returned is link 21, etc. When no active link exists as a south neighbor
                  (in the case of link 19), BAD_INDEX_VALUE is returned. 
                  
    ::

          *------>*------>*------>*------>*
       
       
            
          *--23-->*--24-->*--25-->*--26-->*
    
    
    
          *--19-->*--20-->*--21-->*--22-->*


              
          *------>*------>*------>*------>*
    """

    # First, we find the shape of the horizontal link array given the shape
    # of the raster model grid. In our example, the shape of horizontal links for
    # a grid of 4 rows and 5 columns is 4 rows of horizontal links and 4 columns of
    # horizontal links.
    horizontal_2d_shape = shape_of_horizontal_links(shape)
    
    # Then, we reshape the flattend (1-D) horizontal_link_id array into the shape
    # provided by the shape_of_horizontal_links() function.
    horizontal_2d_array = np.reshape(horizontal_link_ids, horizontal_2d_shape)
    
    # To find south links, we need to shift the IDs in the 2-D array. We first
    # insert a row of BAD_INDEX_VALUE into the top row of the array
    horizontal_link_ids = np.insert(horizontal_2d_array, [0], BAD_INDEX_VALUE, axis=0)
    
    # We find the updated array shape and number of rows for the updated array. 
    row_len = np.shape(horizontal_link_ids)[0]
    
    # To get back to the correct array size (the one found using shape_of_horizontal_links),
    # we delete the last row in the 2-D array
    link_ids  = np.delete(horizontal_link_ids, [row_len-1], axis=0)
    
    # Once we have shifted the 2-D array and removed extra indices, we can flatten
    # the output array to a 1-D array with length of number_of_horizontal_links.
    south_horizontal_neighbors = link_ids.flatten()
    
    return south_horizontal_neighbors   
    
def find_horizontal_west_neighbor(shape, horizontal_link_ids, BAD_INDEX_VALUE=-1):
    """Get IDs of west, horizontal link neighbor

    Parameters
    ----------
    shape : tuple of int
        Shape of grid of nodes.
    horizontal_active_ids : array of int
        Array of all horizontal link ids - MUST BE ARRAY OF LEN(HORIZONTAL_LINKS)
    BAD_INDEX_VALUE: int, optional
        Value assigned to inactive indicies in the array.

    Returns
    -------
    ndarray :
        Link IDs of west horizontal neighbor active links. Length of number_of_horizontal_links.

    Examples
    --------
    >>> from landlab import RasterModelGrid
    >>> from landlab.grid.structured_quad.links import active_link_ids, horizontal_active_link_ids, find_horizontal_west_neighbor

    >>> rmg = RasterModelGrid(4, 5)
    >>> rmg.set_closed_boundaries_at_grid_edges(True, True, True, True)
    >>> status = rmg.node_status 

    >>> active_link_ids = active_link_ids((4,5), status)
    >>> horizontal_active_ids = horizontal_active_link_ids((4,5), active_link_ids)
    >>> find_horizontal_west_neighbor((4,5), horizontal_active_ids)
    array([-1, -1, -1, -1, -1, -1, 20, 21, -1, -1, 24, 25, -1, -1, -1, -1])

    
    Example grid: Indicies are given for active horizontal links in the 4x5 grid space.
                  Only horizontal links are shown, asterisks (*) represent nodes.
                  In this example, link 21 has one west neighbor, link 20. In the 
                  find_horizontal_west_neighbor array, the index of 20 (the west link)
                  corresponds with link 21. Similarly, for link 25, the west neighbor
                  returned is link 24, etc. When no active link exists as a west neighbor
                  (in the case of link 20), BAD_INDEX_VALUE is returned. 

    ::

          *------>*------>*------>*------>*
       
       
            
          *------>*--24-->*--25--->*----->*
    
    
    
          *------>*--20-->*--21-->*------>*


              
          *------>*------>*------>*------>*
    """
    
    # First, we find the shape of the horizontal link array given the shape
    # of the raster model grid. In our example, the shape of horizontal links for
    # a grid of 4 rows and 5 columns is 4 rows of horizontal links and 4 columns of
    # horizontal links. 
    horizontal_2d_shape = shape_of_horizontal_links(shape)
    
    # Then, we reshape the flattend (1-D) horizontal_link_id array into the shape
    # provided by the shape_of_horizontal_links() function. 
    horizontal_2d_array = np.reshape(horizontal_link_ids, horizontal_2d_shape)
     
    # To find west links, we need to shift the IDs in the 2-D array. We insert
    # a column of BAD_INDEX_VALUE into the first column of the array. 
    horizontal_link_ids = np.insert(horizontal_2d_array, [0], BAD_INDEX_VALUE, axis=1)
    
    # We find the updated array shape and number of columns for the updated array. 
    row_len = np.shape(horizontal_link_ids)[1]
    
    # To get back to the correct array size (the one found using shape_of_horizontal_links),
    # we delete the very LAST column of the 2-D array. (Any link final column in the 2-D array
    # cannot be a western neighbor anyway). 
    horizontal_link_ids = np.delete(horizontal_link_ids, [row_len-1], axis=1)
    
    # Once we have shifted the 2-D array and removed extra indices, we can flatten
    # the output array to a 1-D array with length of number_of_horizontal_links.
    west_horizontal_neighbors = horizontal_link_ids.flatten()
    
    return west_horizontal_neighbors


def find_horizontal_north_neighbor(shape, horizontal_link_ids, BAD_INDEX_VALUE=-1):
    """Get IDs of NORTH, horizontal link neighbor

    Parameters
    ----------
    shape : tuple of int
        Shape of grid of nodes.
    horizontal_active_ids : array of int
        Array of all horizontal link ids - MUST BE ARRAY OF LEN(HORIZONTAL_LINKS)
    BAD_INDEX_VALUE: int, optional
        Value assigned to inactive indicies in the array.

    Returns
    -------
    ndarray :
        Link IDs of north horizontal neighbor links. Length of number_of_horizontal_links.

    Examples
    --------
    >>> from landlab import RasterModelGrid
    >>> from landlab.grid.structured_quad.links import active_link_ids, horizontal_active_link_ids, find_horizontal_north_neighbor

    >>> rmg = RasterModelGrid(4, 5)

    >>> active_link_ids = active_link_ids((4,5), rmg.node_status)
    >>> horizontal_active_ids = horizontal_active_link_ids((4,5), active_link_ids)
    >>> find_horizontal_north_neighbor((4,5), horizontal_active_ids)
    array([19, 20, 21, 22, 23, 24, 25, 26, -1, -1, -1, -1, -1, -1, -1, -1])

    Example grid: Indicies are given for active horizontal links in the 4x5 grid space.
                  Only horizontal links are shown, asterisks (*) represent nodes.
                  In this example, link 20 has one north neighbor, link 24. In the 
                  find_horizontal_north_neighbor array, the index of 24 (the north link)
                  corresponds with link 20. Similarly, for link 21, the north neighbor
                  returned is link 25, etc. When no active link exists as a north neighbor
                  (in the case of link 23), BAD_INDEX_VALUE is returned. 

    ::

          *------>*------>*------>*------>*
       
       
            
          *--23-->*--24-->*--25-->*--26-->*
    
    
    
          *--19-->*--20-->*--21-->*--22-->*


              
          *------>*------>*------>*------>*
    """

    # First, we find the shape of the horizontal link array given the shape
    # of the raster model grid. In our example, the shape of horizontal links for
    # a grid of 4 rows and 5 columns is 4 rows of horizontal links and 4 columns of
    # horizontal links.
    horizontal_2d_shape = shape_of_horizontal_links(shape)
    
    # Then, we reshape the flattend (1-D) horizontal_link_id array into the shape
    # provided by the shape_of_horizontal_links() function.
    horizontal_2d_array = np.reshape(horizontal_link_ids, horizontal_2d_shape)
    
    # To find north links, we need to shift the IDs in the 2-D array. We first
    # delete the top row of the array
    horizontal_link_ids = np.delete(horizontal_2d_array, [0], axis=0)
    
    # We find the updated array shape and number of rows for the updated array. 
    row_len = np.shape(horizontal_link_ids)[0]
    
    # To get back to the correct array size (the one found using shape_of_horizontal_links),
    # we insert a row (populated with BAD_INDEX_VALUE_ into the end of the 2-D array.
    link_ids  = np.insert(horizontal_link_ids, [row_len], BAD_INDEX_VALUE, axis=0)
    
    # Once we have shifted the 2-D array and removed extra indices, we can flatten
    # the output array to a 1-D array with length of number_of_horizontal_links.
    north_horizontal_neighbors = link_ids.flatten()
    
    return north_horizontal_neighbors   


def find_horizontal_east_neighbor(shape, horizontal_link_ids, BAD_INDEX_VALUE=-1):
    """Get IDs of east, horizontal link neighbor

    Parameters
    ----------
    shape : tuple of int
        Shape of grid of nodes.
    horizontal_active_ids : array of int
        Array of all horizontal link ids - MUST BE ARRAY OF LEN(HORIZONTAL_LINKS)
    BAD_INDEX_VALUE: int, optional
        Value assigned to inactive indicies in the array.

    Returns
    -------
    ndarray :
        Link IDs of east horizontal neighbor active links. Length of number_of_horizontal_links.

    Examples
    --------
    >>> from landlab import RasterModelGrid
    >>> from landlab.grid.structured_quad.links import active_link_ids, horizontal_active_link_ids, find_horizontal_east_neighbor
    >>>
    >>> rmg = RasterModelGrid(4, 5)
    >>> rmg.set_closed_boundaries_at_grid_edges(True, True, True, True)
    >>> status = rmg.node_status 
    >>> 
    >>> active_link_ids = active_link_ids((4,5), status)
    >>> horizontal_active_ids = horizontal_active_link_ids((4,5), active_link_ids)
    >>> find_horizontal_east_neighbor((4,5), horizontal_active_ids)
    array([-1, -1, -1, -1, 20, 21, -1, -1, 24, 25, -1, -1, -1, -1, -1, -1])

    
    Example grid: Indicies are given for active horizontal links in the 4x5 grid space.
                  Only horizontal links are shown, asterisks (*) represent nodes.
                  In this example, link 20 has one east neighbor, link 21. In the 
                  find_horizontal_east_neighbor array, the index of 21 (the east link)
                  corresponds with link 20. Similarly, for link 24, the east neighbor
                  returned is link 25, etc. When no active link exists as a east neighbor
                  (in the case of link 21), BAD_INDEX_VALUE is returned. 

    ::

          *------>*------>*------>*------>*
       
       
            
          *------>*--24-->*--25--->*----->*
    
    
    
          *------>*--20-->*--21-->*------>*


              
          *------>*------>*------>*------>*
    """
    
    # First, we find the shape of the horizontal link array given the shape
    # of the raster model grid. In our example, the shape of horizontal links for
    # a grid of 4 rows and 5 columns is 4 rows of horizontal links and 4 columns of
    # horizontal links.
    horizontal_2d_shape = shape_of_horizontal_links(shape)
    
    # Then, we reshape the flattend (1-D) horizontal_link_id array into the shape
    # provided by the shape_of_horizontal_links() function.
    horizontal_2d_array = np.reshape(horizontal_link_ids, horizontal_2d_shape)
    
    
    # To find west links, we need to shift the IDs in the 2-D array. We first
    # delete the first column of the array (these values can never be east neighbors
    # anyway.) 
    horizontal_link_ids = np.delete(horizontal_2d_array, [0], axis=1)
    
    # We find the updated array shape and number of columns for the updated array. 
    row_len = np.shape(horizontal_link_ids)[1]
    
    # To get back to the correct array size (the one found using shape_of_horizontal_links),
    # we insert a column of BAD_INDEX_VALUE into the last column spot in the 2-D array.
    link_ids  = np.insert(horizontal_link_ids, [row_len], BAD_INDEX_VALUE, axis=1)
    
    # Once we have shifted the 2-D array and removed extra indices, we can flatten
    # the output array to a 1-D array with length of number_of_horizontal_links.
    east_horizontal_neighbors = link_ids.flatten()
    
    return east_horizontal_neighbors   
               
               
def find_d4_horizontal_neighbors(shape, horizontal_ids, BAD_INDEX_VALUE=-1):
    """Give IDs of all 4 horizontal link neighbors. [S,W,N,E]

    Parameters
    ----------
    shape : tuple of int
        Shape of grid of nodes.
    horizontal_active_ids : array of int
        Array of all horizontal link ids - MUST BE ARRAY OF LEN(HORIZONTAL_LINKS)
    BAD_INDEX_VALUE: int, optional
        Value assigned to inactive indicies in the array.

    Returns
    -------
    ndarray :
        Array of 4 horizontal link neighbors for a given link ID. Returned in
        [S, W, N, E]. 

    Examples
    --------
    >>> from landlab import RasterModelGrid
    >>> from landlab.grid.structured_quad.links import active_link_ids, horizontal_active_link_ids, find_d4_horizontal_neighbors
    >>> rmg = RasterModelGrid(4, 5)
    >>> active_link_ids = active_link_ids((4,5), rmg.node_status)
    >>> horizontal_active_ids = horizontal_active_link_ids((4,5), active_link_ids)
    >>> find_d4_horizontal_neighbors((4,5), horizontal_active_ids)
    array([[-1, -1, 19, -1],
           [-1, -1, 20, -1],
           [-1, -1, 21, -1],
           [-1, -1, 22, -1],
           [-1, -1, 23, 20],
           [-1, 19, 24, 21],
           [-1, 20, 25, 22],
           [-1, 21, 26, -1],
           [19, -1, -1, 24],
           [20, 23, -1, 25],
           [21, 24, -1, 26],
           [22, 25, -1, -1],
           [23, -1, -1, -1],
           [24, -1, -1, -1],
           [25, -1, -1, -1],
           [26, -1, -1, -1]])

    
    Example grid: Indicies are given for active horizontal links in the 4x5 grid space.
                  Only horizontal links are shown, asterisks (*) represent nodes.
                  
                  In this example, link 20 has three neighbors. (links 19, 21, and 24).
                  This function looks for S, W, N, E neighbors and returns an array of
                  indices for each link. For link 20, the returned array would be
                  [-1, 19, 24, 21]. 

    ::

          *------>*------>*------>*------>*
       
       
            
          *------>*--24-->*------>*------>*
    
    
    
          *--19-->*--20-->*--21-->*------>*


              
          *------>*------>*------>*------>*
    """
    ### First we find *south* neighbors...
    south = find_horizontal_south_neighbor(shape, horizontal_ids, BAD_INDEX_VALUE)
    
    ### Then *west* neighbors...
    west = find_horizontal_west_neighbor(shape, horizontal_ids, BAD_INDEX_VALUE) 

    ### Then *north* neighbors...
    north = find_horizontal_north_neighbor(shape, horizontal_ids, BAD_INDEX_VALUE)

    ### Finally, *east* neighbors...
    east = find_horizontal_east_neighbor(shape, horizontal_ids, BAD_INDEX_VALUE)
   
    ### Combine all 4 neighbor arrays into one large array (4 x len_horizontal_links)
    neighbor_array = np.array([south, west, north, east])
    
    ### Transpose the 4 neighbor arrays into a (len_horizontal_links x 4) array.
    neighbor_array = np.transpose(neighbor_array)
    
    ### Output neighbor array. For each input ID, returns [S,W,N,E]
    return neighbor_array            
       
       
def find_d4_horizontal_neighbors_active(shape, horizontal_ids, BAD_INDEX_VALUE=-1):
    """Give IDs of all 4 horizontal link neighbors. [S,W,N,E]

    Parameters
    ----------
    shape : tuple of int
        Shape of grid of nodes.
    horizontal_active_ids : array of int
        Array of all horizontal link ids - MUST BE ARRAY OF LEN(HORIZONTAL_LINKS)
    BAD_INDEX_VALUE: int, optional
        Value assigned to inactive indicies in the array.

    Returns
    -------
    ndarray :
        Array of 4 horizontal link neighbors for a given link ID. Returned in
        [S, W, N, E]. Returns array for only ACTIVE horizontal links.

    Examples
    --------
    >>> from landlab import RasterModelGrid
    >>> from landlab.grid.structured_quad.links import active_link_ids, horizontal_active_link_ids, find_d4_horizontal_neighbors_active

    >>> rmg = RasterModelGrid(4, 5)

    >>> active_link_ids = active_link_ids((4,5), rmg.node_status)
    >>> horizontal_active_ids = horizontal_active_link_ids((4,5), active_link_ids)
    >>> find_d4_horizontal_neighbors_active((4,5), horizontal_active_ids)
    array([[-1, -1, 23, 20],
           [-1, 19, 24, 21],
           [-1, 20, 25, 22],
           [-1, 21, 26, -1],
           [19, -1, -1, 24],
           [20, 23, -1, 25],
           [21, 24, -1, 26],
           [22, 25, -1, -1]])


    Example grid: Indicies are given for active horizontal links in the 4x5 grid space.
                  Only horizontal links are shown, asterisks (*) represent nodes.
                  
                  In this example, link 20 has three neighbors. (links 19, 21, and 24).
                  This function looks for S, W, N, E neighbors and returns an array of
                  indices for each link. For link 20, the returned array would be
                  [-1, 19, 24, 21]. 

    ::

          *------>*------>*------>*------>*
       
       
            
          *------>*--24-->*------>*------>*
    
    
    
          *--19-->*--20-->*--21-->*------>*


              
          *------>*------>*------>*------>*
    """
    # To do this we simply call the find_d4_horizontal_neighbors() function 
    # which gives the neighbors for ALL horizontal links in an array, even 
    # inactive links. 
    d4_neigh = find_d4_horizontal_neighbors(shape, horizontal_ids,
                                            BAD_INDEX_VALUE)
    
    # Now we will just focus on indices that are ACTIVE...
    active_links = np.where(horizontal_ids != BAD_INDEX_VALUE)
    
    # Clip our initial array into a smaller one with just active neighbors
    neighbor_array = d4_neigh[active_links]
    
    # Output neighbor array. For each input ID, returns [S,W,N,E]
    return neighbor_array                                  


def find_vertical_south_neighbor(shape, vertical_link_ids, BAD_INDEX_VALUE=-1):
    """Link IDs of south, vertical link neighbor

    Parameters
    ----------
    shape : tuple of int
        Shape of grid of nodes.
    vertical_link_ids : array of int
        Array of all vertical link ids - MUST BE ARRAY OF LEN(VERTICAL_LINKS)
    BAD_INDEX_VALUE: int, optional
        Value assigned to inactive indicies in the array.

    Returns
    -------
    ndarray :
        Link IDs of *south* vertical neighbor active links. Length of
        number_of_vertical_links.

    Examples
    --------
    >>> from landlab import RasterModelGrid
    >>> from landlab.grid.structured_quad.links import active_link_ids, vertical_active_link_ids, find_vertical_south_neighbor

    >>> rmg = RasterModelGrid(4, 5)

    >>> active_link_ids = active_link_ids((4,5), rmg.node_status)
    >>> vertical_active_ids = vertical_active_link_ids((4,5), active_link_ids)
    >>> find_vertical_south_neighbor((4,5), vertical_active_ids)
    array([-1, -1, -1, -1, -1, -1,  1,  2,  3, -1, -1,  6,  7,  8, -1])

    Example grid: Indicies are given for active vertical links in the 4x5 grid space.
                  Only vertical links are shown, asterisks (*) represent nodes.
                  In this example, link 6 has one south neighbor, link 1. In the 
                  find_vertical_north_neighbor array, the index of 1 (the south link)
                  corresponds with link 6. Similarly, for link 13, the south neighbor
                  returned is link 8, etc. When no active link exists as a east neighbor
                  (in the case of link 2), BAD_INDEX_VALUE is returned. 
    ::

          *       *       *       *       *
          ^       ^       ^       ^       ^   
          |       11      12      13      |      
          |       |       |       |       |         
          *       *       *       *       *
          ^       ^       ^       ^       ^       
          |       6       7       8       |       
          |       |       |       |       |       
          *       *       *       *       *
          ^       ^       ^       ^       ^   
          |       1       2       3       |       
          |       |       |       |       |       
          *       *       *       *       *
    """
    # First, we find the shape of the vertical link array given the shape
    # of the raster model grid. In our example, the shape of vertical links for
    # a grid of 4 rows and 5 columns is 3 rows of vertical links and 5 columns of
    # vertical links.   
    vertical_2d_shape = shape_of_vertical_links(shape)
    
    # Then, we reshape the flattend (1-D) vertical_link_id array into the shape
    # provided by the shape_of_vertical_links() function.
    vertical_2d_array = np.reshape(vertical_link_ids, vertical_2d_shape)

    # To find south links, we need to shift the IDs in the 2-D array. We insert
    # a row of BAD_INDEX_VALUE into the top row of the 2-D array
    link_ids  = np.insert(vertical_2d_array, [0],  BAD_INDEX_VALUE, axis=0)
    
    # We find the updated array shape and number of rows for the updated array. 
    row_len = np.shape(link_ids)[0]
    
    # To get back to the correct array size (the one found using shape_of_vertical_links),
    # we delete a the last row of the 2-D array. 
    vertical_link_ids = np.delete(link_ids, [row_len-1], axis=0)

    # Once we have shifted the 2-D array and removed extra indices, we can flatten
    # the output array to a 1-D array with length of number_of_vertical_links.
    south_vertical_neighbors = vertical_link_ids.flatten()
    
    return south_vertical_neighbors      

def find_vertical_west_neighbor(shape, vertical_link_ids, BAD_INDEX_VALUE=-1):
    """Link IDs of west, vertical link neighbor

    Parameters
    ----------
    shape : tuple of int
        Shape of grid of nodes.
    vertical_link_ids : array of int
        Array of all vertical link ids- MUST BE ARRAY OF LEN(VERTICAL_LINKS)
    BAD_INDEX_VALUE: int, optional
        Value assigned to inactive indicies in the array.

    Returns
    -------
    ndarray :
        Link IDs of *west* vertical neighbor active links. Length of number_of_vertical_links.

    Examples
    --------
    >>> from landlab import RasterModelGrid
    >>> from landlab.grid.structured_quad.links import active_link_ids, vertical_active_link_ids, find_vertical_west_neighbor

    >>> rmg = RasterModelGrid(4, 5)
    >>> rmg.set_closed_boundaries_at_grid_edges(True, True, True, True)
    >>> status = rmg.node_status 

    >>> active_link_ids = active_link_ids((4,5), status)
    >>> vertical_active_ids = vertical_active_link_ids((4,5), active_link_ids)
    >>> find_vertical_west_neighbor((4,5), vertical_active_ids)
    array([-1, -1, -1, -1, -1, -1, -1,  6,  7,  8, -1, -1, -1, -1, -1])


    Example grid: Indicies are given for active vertical links in the 4x5 grid space.
                  Only vertical links are shown, asterisks (*) represent nodes.
                  In this example, link 7 has one west neighbor, link 6. In the 
                  find_vertical_west_neighbor array, the index of 6 (the west link)
                  corresponds with link 7. Similarly, for link 8, the west neighbor
                  returned is link 7, etc. When no active link exists as a west neighbor
                  (in the case of link 6), BAD_INDEX_VALUE is returned. 

    ::

          *       *       *       *       *
          ^       ^       ^       ^       ^   
          |       |       |       |       |      
          |       |       |       |       |         
          *       *       *       *       *
          ^       ^       ^       ^       ^       
          |       6       7       8       |       
          |       |       |       |       |       
          *       *       *       *       *
          ^       ^       ^       ^       ^   
          |       |       |       |       |       
          |       |       |       |       |       
          *       *       *       *       *
    """
    # First, we find the shape of the vertical link array given the shape
    # of the raster model grid. In our example, the shape of vertical links for
    # a grid of 4 rows and 5 columns is 3 rows of vertical links and 5 columns of
    # vertical links. 
    vertical_2d_shape = shape_of_vertical_links(shape)
    
    # Then, we reshape the flattend (1-D) vertical_link_id array into the shape
    # provided by the shape_of_vertical_links() function. 
    vertical_2d_array = np.reshape(vertical_link_ids, vertical_2d_shape)
    
    # To find west links, we need to shift the IDs in the 2-D array. We insert
    # a column of BAD_INDEX_VALUE into the first column of the array. 
    vertical_link_ids = np.insert(vertical_2d_array, [0], BAD_INDEX_VALUE, axis=1)
    
    # We find the updated array shape and number of columns for the updated array. 
    row_len = np.shape(vertical_link_ids)[1]
    
    # To get back to the correct array size (the one found using shape_of_vertical_links),
    # we delete the very LAST column of the 2-D array. (Any link final column in the 2-D array
    # cannot be a western neighbor anyway). 
    vertical_link_ids = np.delete(vertical_link_ids, [row_len-1], axis=1)
    
    # Once we have shifted the 2-D array and removed extra indices, we can flatten
    # the output array to a 1-D array with length of number_of_vertical_links.
    west_vertical_neighbors = vertical_link_ids.flatten()
    
    return west_vertical_neighbors


def find_vertical_north_neighbor(shape, vertical_link_ids, BAD_INDEX_VALUE=-1):
    """Link IDs of north, vertical link neighbor

    Parameters
    ----------
    shape : tuple of int
        Shape of grid of nodes.
    vertical_link_ids : array of int
        Array of all vertical link ids- MUST BE ARRAY OF LEN(VERTICAL_LINKS)
    BAD_INDEX_VALUE: int, optional
        Value assigned to inactive indicies in the array.

    Returns
    -------
    ndarray :
        Link IDs of *north* vertical neighbor active links. Length of number_of_vertical_links.

    Examples
    --------
    >>> from landlab import RasterModelGrid
    >>> from landlab.grid.structured_quad.links import active_link_ids, vertical_active_link_ids, find_vertical_north_neighbor

    >>> rmg = RasterModelGrid(4, 5)

    >>> active_link_ids = active_link_ids((4,5), rmg.node_status)
    >>> vertical_active_ids = vertical_active_link_ids((4,5), active_link_ids)
    >>> find_vertical_north_neighbor((4,5), vertical_active_ids)
    array([-1,  6,  7,  8, -1, -1, 11, 12, 13, -1, -1, -1, -1, -1, -1])


    Example grid: Indicies are given for active vertical links in the 4x5 grid space.
                  Only vertical links are shown, asterisks (*) represent nodes.
                  In this example, link 1 has one north neighbor, link 6. In the 
                  find_vertical_north_neighbor array, the index of 6 (the north link)
                  corresponds with link 1. Similarly, for link 8, the north neighbor
                  returned is link 13, etc. When no active link exists as a east neighbor
                  (in the case of link 11), BAD_INDEX_VALUE is returned. 

    ::

          *       *       *       *       *
          ^       ^       ^       ^       ^   
          |       11      12      13      |      
          |       |       |       |       |         
          *       *       *       *       *
          ^       ^       ^       ^       ^       
          |       6       7       8       |       
          |       |       |       |       |       
          *       *       *       *       *
          ^       ^       ^       ^       ^   
          |       1       2       3       |       
          |       |       |       |       |       
          *       *       *       *       *
    """
    # First, we find the shape of the vertical link array given the shape
    # of the raster model grid. In our example, the shape of vertical links for
    # a grid of 4 rows and 5 columns is 3 rows of vertical links and 5 columns of
    # vertical links.   
    vertical_2d_shape = shape_of_vertical_links(shape)
    
    # Then, we reshape the flattend (1-D) vertical_link_id array into the shape
    # provided by the shape_of_vertical_links() function.
    vertical_2d_array = np.reshape(vertical_link_ids, vertical_2d_shape)

    # To find north links, we need to shift the IDs in the 2-D array. We first
    # delete the first row of the array.
    vertical_link_ids = np.delete(vertical_2d_array, [0], axis=0)

    # We find the updated array shape and number of rows for the updated array. 
    row_len = np.shape(vertical_link_ids)[0]
    
    # To get back to the correct array size (the one found using shape_of_vertical_links),
    # we insert a row (populated with BAD_INDEX_VALUE) into the end of the 2-D array. 
    link_ids  = np.insert(vertical_link_ids, [row_len],  BAD_INDEX_VALUE, axis=0)
    
    # Once we have shifted the 2-D array and removed extra indices, we can flatten
    # the output array to a 1-D array with length of number_of_vertical_links.
    north_vertical_neighbors = link_ids.flatten()
    
    return north_vertical_neighbors
    
    
def find_vertical_east_neighbor(shape, vertical_link_ids, BAD_INDEX_VALUE=-1):
    """Link IDs of east, vertical link neighbor

    Parameters
    ----------
    shape : tuple of int
        Shape of grid of nodes.
    vertical_link_ids : array of int
        Array of all vertical link ids - MUST BE ARRAY OF LEN(VERTICAL_LINKS)
    BAD_INDEX_VALUE: int, optional
        Value assigned to inactive indicies in the array.

    Returns
    -------
    ndarray :
        Link IDs of *east* vertical neighbor active links. Length of number_of_vertical_links.

    Examples
    --------
    >>> from landlab import RasterModelGrid
    >>> from landlab.grid.structured_quad.links import active_link_ids, vertical_active_link_ids, find_vertical_east_neighbor

    >>> rmg = RasterModelGrid(4, 5)
    >>> rmg.set_closed_boundaries_at_grid_edges(True, True, True, True)
    >>> status = rmg.node_status 

    >>> active_link_ids = active_link_ids((4,5), status)
    >>> vertical_active_ids = vertical_active_link_ids((4,5), active_link_ids)
    >>> find_vertical_east_neighbor((4,5), vertical_active_ids)
    array([-1, -1, -1, -1, -1,  6,  7,  8, -1, -1, -1, -1, -1, -1, -1])


    Example grid: Indicies are given for active vertical links in the 4x5 grid space.
                  Only vertical links are shown, asterisks (*) represent nodes.
                  In this example, link 6 has one east neighbor, link 7. In the 
                  find_vertical_east_neighbor array, the index of 7 (the east link)
                  corresponds with link 6. Similarly, for link 7, the east neighbor
                  returned is link 8, etc. When no active link exists as a east neighbor
                  (in the case of link 8), BAD_INDEX_VALUE is returned. 
                
    ::

          *       *       *       *       *
          ^       ^       ^       ^       ^   
          |       |       |       |       |      
          |       |       |       |       |         
          *       *       *       *       *
          ^       ^       ^       ^       ^       
          |       6       7       8       |       
          |       |       |       |       |       
          *       *       *       *       *
          ^       ^       ^       ^       ^   
          |       |       |       |       |       
          |       |       |       |       |       
          *       *       *       *       *
    """
    # First, we find the shape of the vertical link array given the shape
    # of the raster model grid. In our example, the shape of vertical links for
    # a grid of 4 rows and 5 columns is 3 rows of vertical links and 5 columns of
    # vertical links. 
    vertical_2d_shape = shape_of_vertical_links(shape)
    
    # Then, we reshape the flattend (1-D) vertical_link_id array into the shape
    # provided by the shape_of_vertical_links() function.
    vertical_2d_array = np.reshape(vertical_link_ids, vertical_2d_shape)
    
    
    # To find east links, we need to shift the IDs in the 2-D array. We first
    # delete the first column of the array.
    vertical_link_ids = np.delete(vertical_2d_array, [0], axis=1)
    
    # We find the updated array shape and number of columns for the updated array. 
    row_len = np.shape(vertical_link_ids)[1]
    
    # To get back to the correct array size (the one found using shape_of_vertical_links),
    # we insert a column (populated with BAD_INDEX_VALUE) into the end of the 2-D array. 
    link_ids  = np.insert(vertical_link_ids, [row_len],  BAD_INDEX_VALUE, axis=1)
    
    # Once we have shifted the 2-D array and removed extra indices, we can flatten
    # the output array to a 1-D array with length of number_of_vertical_links.
    east_vertical_neighbors = link_ids.flatten()
    
    return east_vertical_neighbors    
    
    
def find_d4_vertical_neighbors(shape, vertical_ids, BAD_INDEX_VALUE=-1):
    """Give IDs of all 4 vertical link neighbors. [S,W,N,E]

    Parameters
    ----------
    shape : tuple of int
        Shape of grid of nodes.
    vertical_active_ids : array of int
        Array of all vertical link ids - MUST BE ARRAY OF LEN(VERTICAL_LINKS)
    BAD_INDEX_VALUE: int, optional
        Value assigned to inactive indicies in the array.

    Returns
    -------
    ndarray :
        Array of 4 vertical link neighbors for a given link ID. Returned in
        [S, W, N, E]. 

    Examples
    --------
    >>> from landlab import RasterModelGrid
    >>> from landlab.grid.structured_quad.links import active_link_ids, vertical_active_link_ids, find_d4_vertical_neighbors

    >>> rmg = RasterModelGrid(4, 5)

    >>> active_link_ids = active_link_ids((4,5), rmg.node_status)
    >>> vertical_active_ids = vertical_active_link_ids((4,5), active_link_ids)
    >>> find_d4_vertical_neighbors((4,5), vertical_active_ids)
    array([[-1, -1, -1,  1],
           [-1, -1,  6,  2],
           [-1,  1,  7,  3],
           [-1,  2,  8, -1],
           [-1,  3, -1, -1],
           [-1, -1, -1,  6],
           [ 1, -1, 11,  7],
           [ 2,  6, 12,  8],
           [ 3,  7, 13, -1],
           [-1,  8, -1, -1],
           [-1, -1, -1, 11],
           [ 6, -1, -1, 12],
           [ 7, 11, -1, 13],
           [ 8, 12, -1, -1],
           [-1, 13, -1, -1]])
    
    
    Example grid: Indicies are given for active vertical links in the 4x5 grid space.
                  Only vertical links are shown, asterisks (*) represent nodes.
                  
                  In this example, link 7 has four neighbors. (links 2, 6, 12 and 8 in 
                  S, W, N, E order)
                  
                  This function looks for S, W, N, E neighbors and returns an array of
                  indices for each link. For link 8, the returned array would be
                  [3, 7, 13, -1]. 

    ::

          *       *       *       *       *
          ^       ^       ^       ^       ^   
          |       11      12     13       |      
          |       |       |       |       |         
          *       *       *       *       *
          ^       ^       ^       ^       ^       
          |       6       7       8       |       
          |       |       |       |       |       
          *       *       *       *       *
          ^       ^       ^       ^       ^   
          |       1       2       3       |       
          |       |       |       |       |       
          *       *       *       *       *
    """
    south = find_vertical_south_neighbor(shape, vertical_ids, BAD_INDEX_VALUE)
    west = find_vertical_west_neighbor(shape, vertical_ids, BAD_INDEX_VALUE) 
    north = find_vertical_north_neighbor(shape, vertical_ids, BAD_INDEX_VALUE)
    east = find_vertical_east_neighbor(shape, vertical_ids, BAD_INDEX_VALUE)   
    neighbor_array = np.array([south, west, north, east])
    neighbor_array = np.transpose(neighbor_array)
    return neighbor_array            


def find_d4_vertical_neighbors_active(shape, vertical_ids, BAD_INDEX_VALUE=-1):
    """Give IDs of all 4 vertical link neighbors. [S,W,N,E]

    Parameters
    ----------
    shape : tuple of int
        Shape of grid of nodes.
    vertical_active_ids : array of int
        Array of all vertical link ids - MUST BE ARRAY OF LEN(VERTICAL_LINKS)
    BAD_INDEX_VALUE: int, optional
        Value assigned to inactive indicies in the array.

    Returns
    -------
    ndarray :
        Array of 4 vertical link neighbors for a given ACTIVE link ID. Returned in
        [S, W, N, E]. 

    Examples
    --------
    >>> from landlab import RasterModelGrid
    >>> from landlab.grid.structured_quad.links import active_link_ids, vertical_active_link_ids, find_d4_vertical_neighbors_active

    >>> rmg = RasterModelGrid(4, 5)

    >>> active_link_ids = active_link_ids((4,5), rmg.node_status)
    >>> vertical_active_ids = vertical_active_link_ids((4,5), active_link_ids)
    >>> find_d4_vertical_neighbors_active((4,5), vertical_active_ids)
    array([[-1, -1,  6,  2],
           [-1,  1,  7,  3],
           [-1,  2,  8, -1],
           [ 1, -1, 11,  7],
           [ 2,  6, 12,  8],
           [ 3,  7, 13, -1],
           [ 6, -1, -1, 12],
           [ 7, 11, -1, 13],
           [ 8, 12, -1, -1]])

    
    Example grid: Indicies are given for active horizontal links in the 4x5 grid space.
                  Only horizontal links are shown, asterisks (*) represent nodes.
                  
                  In this example, link 7 has four neighbors. (links 2, 6, 12 and 8 in 
                  S, W, N, E order)
                  
                  This function looks for S, W, N, E neighbors and returns an array of
                  indices for each link. For link 8, the returned array would be
                  [3, 7, 13, -1]. 

    ::

          *       *       *       *       *
          ^       ^       ^       ^       ^   
          |       11      12     13       |      
          |       |       |       |       |         
          *       *       *       *       *
          ^       ^       ^       ^       ^       
          |       6       7       8       |       
          |       |       |       |       |       
          *       *       *       *       *
          ^       ^       ^       ^       ^   
          |       1       2       3       |       
          |       |       |       |       |       
          *       *       *       *       *
    """
    # To do this we simply call the find_d4_vertical_neighbors() function 
    # which gives the neighbors for ALL vertical links in an array, even 
    # inactive links. 
    d4_all_neighbors = find_d4_vertical_neighbors(shape, vertical_ids, BAD_INDEX_VALUE)

    # Now we will just focus on indices that are ACTIVE...
    active_links = np.where(vertical_ids != BAD_INDEX_VALUE)
    
    # Clip our initial array into a smaller one with just active neighbors
    neighbor_array = d4_all_neighbors[active_links]
    
    # Output neighbor array. For each input ID, returns [S,W,N,E]
    return neighbor_array                

def bottom_edge_horizontal_ids(shape):
    """Link IDs of bottom edge horizontal links

    Parameters
    ----------
    shape : tuple of int
        Shape of grid, given as (rows, columns) of nodes.

    Returns
    -------
    ndarray :
        Link IDs of bottom edge horizontal links. Length is (rmg.number_of_columns-1)

    Examples
    --------
    >>> from landlab import RasterModelGrid
    >>> from landlab.grid.structured_quad.links import bottom_edge_horizontal_ids

    >>> rmg = RasterModelGrid(4, 5)
    >>> shape = rmg.shape

    >>> bottom_edge_horizontal_ids(shape)
    array([15, 16, 17, 18])


    Example grid: Indicies are given for horizontal links along the bottom edge of
                  the 4x5 grid space. Only horizontal links are shown, asterisks (*) 
                  represent nodes.

    ::

          *------>*------>*------>*------>*
       
       
            
          *------>*------>*------>*------>*
    
    
    
          *------>*------>*------>*------>*


              
          *--15-->*--16-->*--17-->*--18-->*
    """
    
    #First, we find all horizontal link ids for the RasterModelGrid shape.
    horizontal_id_array = horizontal_link_ids(shape)    
    
    # Then we slice the first column and return it. This has our bottom edge
    # horizontal ids. This array should be equal in length to (number of columns-1)
    bottom_edge_horizontal_ids = horizontal_id_array[0]

    return bottom_edge_horizontal_ids

def left_edge_horizontal_ids(shape):
    """Link IDs of left edge horizontal links

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
    >>> from landlab import RasterModelGrid
    >>> from landlab.grid.structured_quad.links import left_edge_horizontal_ids

    >>> rmg = RasterModelGrid(4, 5)
    >>> shape = rmg.shape

    >>> left_edge_horizontal_ids(shape)
    array([15, 19, 23, 27])


    Example grid: Indicies are given for horizontal links along the left edge of
                  the 4x5 grid space. Only horizontal links are shown, asterisks (*) 
                  represent nodes.

    ::
          *--27-->*------>*------>*------>*
       
       
            
          *--23-->*------>*------>*------>*
    
    
    
          *--19-->*------>*------>*------>*


              
          *--15-->*------>*------>*------>*
    """
    
    #First, we find all horizontal link ids for the RasterModelGrid shape.
    horizontal_id_array = horizontal_link_ids(shape)    
    
    # Then we slice the first column and return it. This has our left edge
    # horizontal ids. This array should be equal in length to (number of rows)
    left_edge_horizontal_ids = horizontal_id_array[:,0]

    return left_edge_horizontal_ids

def top_edge_horizontal_ids(shape):
    """Link IDs of top edge horizontal links

    Parameters
    ----------
    shape : tuple of int
        Shape of grid, given as (rows, columns) of nodes.

    Returns
    -------
    ndarray :
        Link IDs of top edge horizontal links. Length is (rmg.number_of_columns-1)

    Examples
    --------
    >>> from landlab import RasterModelGrid
    >>> from landlab.grid.structured_quad.links import top_edge_horizontal_ids

    >>> rmg = RasterModelGrid(4, 5)
    >>> shape = rmg.shape

    >>> top_edge_horizontal_ids(shape)
    array([27, 28, 29, 30])


    Example grid: Indicies are given for horizontal links along the top edge of
                  the 4x5 grid space. Only horizontal links are shown, asterisks (*) 
                  represent nodes.

    ::

          *--27-->*--28-->*--29-->*--30-->*
       
       
            
          *------>*------>*------>*------>*
    
    
    
          *------>*------>*------>*------>*


              
          *------>*------>*------>*------>*
    """
    #First, we find all horizontal link ids for the RasterModelGrid shape.
    horizontal_id_array = horizontal_link_ids(shape)    
    
    # Then we slice the first column and return it. This has our top edge
    # horizontal ids. This array should be equal in length to (number of columns-1)
    top_edge_horizontal_ids = horizontal_id_array[(shape[0]-1)]

    return top_edge_horizontal_ids

def right_edge_horizontal_ids(shape):
    """Link IDs of right edge horizontal links

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
    >>> from landlab import RasterModelGrid
    >>> from landlab.grid.structured_quad.links import right_edge_horizontal_ids

    >>> rmg = RasterModelGrid(4, 5)
    >>> shape = rmg.shape

    >>> right_edge_horizontal_ids(shape)
    array([18, 22, 26, 30])


    Example grid: Indicies are given for horizontal links along the right edge of
                  the 4x5 grid space. Only horizontal links are shown, asterisks (*) 
                  represent nodes.

    ::
          *------>*------>*------>*--30-->*
       
       
            
          *------>*------>*------>*--26-->*
    
    
    
          *------>*------>*------>*--22-->*


              
          *------>*------>*------>*--18-->*
    """
    
    #First, we find all horizontal link ids for the RasterModelGrid shape.
    horizontal_id_array = horizontal_link_ids(shape)    
    
    # Then we slice the last column and return it. This has our right edge
    # horizontal ids. This array should be equal in length to (number of rows - 1)
    right_edge_horizontal_ids = horizontal_id_array[:,(shape[0]-1)]

    return right_edge_horizontal_ids
                          
def bottom_edge_vertical_ids(shape):
    """Link IDs of bottom edge vertical links

    Parameters
    ----------
    shape : tuple of int
        Shape of grid, given as (rows, columns) of nodes.

    Returns
    -------
    ndarray :
        Link IDs of bottom edge vertical links. Length is (rmg.number_of_columns)

    Examples
    --------
    >>> from landlab import RasterModelGrid
    >>> from landlab.grid.structured_quad.links import bottom_edge_vertical_ids

    >>> rmg = RasterModelGrid(4, 5)
    >>> shape = rmg.shape

    >>> bottom_edge_vertical_ids(shape)
    array([0, 1, 2, 3, 4])


    Example grid: Indicies are given for vertical links along the bottom edge of
                  the 4x5 grid space. Only vertical links are shown, asterisks (*) 
                  represent nodes.

    ::

          *       *       *       *       *
          ^       ^       ^       ^       ^   
          |       |       |       |       |      
          |       |       |       |       |         
          *       *       *       *       *
          ^       ^       ^       ^       ^       
          |       |       |       |       |       
          |       |       |       |       |       
          *       *       *       *       *
          ^       ^       ^       ^       ^   
          0       1       2       3       4       
          |       |       |       |       |       
          *       *       *       *       *
    """
    
    #First, we find all vertical link ids for the RasterModelGrid shape.
    vertical_id_array = vertical_link_ids(shape)    
    
    # Then we slice the first column and return it. This has our bottom edge
    # vertical ids. This array should be equal in length to (number of columns)
    bottom_edge_vertical_ids = vertical_id_array[0]

    return bottom_edge_vertical_ids

def left_edge_vertical_ids(shape):
    """Link IDs of left edge vertical links

    Parameters
    ----------
    shape : tuple of int
        Shape of grid, given as (rows, columns) of nodes.

    Returns
    -------
    ndarray :
        Link IDs of left edge vertical links. Length is (rmg.number_of_rows - 1)

    Examples
    --------
    >>> from landlab import RasterModelGrid
    >>> from landlab.grid.structured_quad.links import left_edge_vertical_ids

    >>> rmg = RasterModelGrid(4, 5)
    >>> shape = rmg.shape

    >>> left_edge_vertical_ids(shape)
    array([ 0,  5, 10])


    Example grid: Indicies are given for vertical links along the left edge of
                  the 4x5 grid space. Only vertical links are shown, asterisks (*) 
                  represent nodes.

    ::

          *       *       *       *       *
          ^       ^       ^       ^       ^   
          10      |       |       |       |      
          |       |       |       |       |         
          *       *       *       *       *
          ^       ^       ^       ^       ^       
          5       |       |       |       |       
          |       |       |       |       |       
          *       *       *       *       *
          ^       ^       ^       ^       ^   
          0       |       |       |       |       
          |       |       |       |       |       
          *       *       *       *       *
    """
    
    #First, we find all vertical link ids for the RasterModelGrid shape.
    vertical_id_array = vertical_link_ids(shape)    
    
    # Then we slice the first column and return it. This has our left edge
    # vertical ids. This array should be equal in length to (number of rows - 1)
    left_edge_vertical_ids = vertical_id_array[:,0]

    return left_edge_vertical_ids

def top_edge_vertical_ids(shape):
    """Link IDs of top edge vertical links

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
    >>> from landlab import RasterModelGrid
    >>> from landlab.grid.structured_quad.links import top_edge_vertical_ids

    >>> rmg = RasterModelGrid(4, 5)
    >>> shape = rmg.shape

    >>> top_edge_vertical_ids(shape)
    array([10, 11, 12, 13, 14])


    Example grid: Indicies are given for vertical links along the top edge of
                  the 4x5 grid space. Only vertical links are shown, asterisks (*) 
                  represent nodes.

    ::

          *       *       *       *       *
          ^       ^       ^       ^       ^   
          10      11      12      13      14      
          |       |       |       |       |         
          *       *       *       *       *
          ^       ^       ^       ^       ^       
          |       |       |       |       |       
          |       |       |       |       |       
          *       *       *       *       *
          ^       ^       ^       ^       ^   
          |       |       |       |       |      
          |       |       |       |       |       
          *       *       *       *       *
    """
    #First, we find all vertical link ids for the RasterModelGrid shape.
    vertical_id_array = vertical_link_ids(shape)    
    
    # Then we slice the first column and return it. This has our top edge
    # vertical ids. This array should be equal in length to (number of columns)
    top_edge_vertical_ids = vertical_id_array[(shape[0]-2)]

    return top_edge_vertical_ids

def right_edge_vertical_ids(shape):
    """Link IDs of right edge vertical links

    Parameters
    ----------
    shape : tuple of int
        Shape of grid, given as (rows, columns) of nodes.

    Returns
    -------
    ndarray :
        Link IDs of left edge vertical links. Length is (rmg.number_of_rows - 1)

    Examples
    --------
    >>> from landlab import RasterModelGrid
    >>> from landlab.grid.structured_quad.links import right_edge_vertical_ids

    >>> rmg = RasterModelGrid(4, 5)
    >>> shape = rmg.shape

    >>> right_edge_vertical_ids(shape)
    array([ 4,  9, 14])


    Example grid: Indicies are given for vertical links along the right edge of
                  the 4x5 grid space. Only vertical links are shown, asterisks (*) 
                  represent nodes.

    ::

          *       *       *       *       *
          ^       ^       ^       ^       ^   
          |       |       |       |       14      
          |       |       |       |       |         
          *       *       *       *       *
          ^       ^       ^       ^       ^       
          |       |       |       |       9       
          |       |       |       |       |       
          *       *       *       *       *
          ^       ^       ^       ^       ^   
          |       |       |       |       4       
          |       |       |       |       |       
          *       *       *       *       *
    """
    
    #First, we find all vertical link ids for the RasterModelGrid shape.
    vertical_id_array = vertical_link_ids(shape)    
    
    # Then we slice the last column and return it. This has our right edge
    # vertical ids. This array should be equal in length to (number of rows - 1)
    right_edge_vertical_ids = vertical_id_array[:,(shape[1]-1)]

    return right_edge_vertical_ids


class StructuredQuadLinkGrid(LinkGrid):
    def __init__(self, shape):
        link_ends = (node_id_at_link_start(shape), node_id_at_link_end(shape))
        number_of_nodes = np.prod(shape)
        LinkGrid.__init__(self, link_ends, number_of_nodes)

#! /usr/bin/env python


import numpy as np


_BAD_INDEX_VALUE = np.iinfo(np.int).max


def node_count(shape):
    assert(len(shape) == 2)
    return shape[0] * shape[1]


def cell_count(shape):
    assert(len(shape) == 2)
    return (shape[0] - 2) * (shape[1] - 2)


def active_cell_count(shape):
    return cell_count(shape)


def active_link_count(shape):
    assert(len(shape) == 2)
    return link_count(shape) - 2 * (shape[0] - 1) + 2 * (shape[1] - 1)


def link_count(shape):
    """
    Returns the number of links in a structured grid with dimensions, *shape*.
    This is the number of to-links and from-links, not the total of the two.

    >>> link_count((3,2))
    7
    """
    assert(len(shape) == 2)
    return shape[1] * (shape[0] - 1) + shape[0] * (shape[1] - 1)


def boundary_cell_count(shape):
    assert(len(shape) == 2)
    return 2 * (shape[0] - 2) + 2 * (shape[1] - 2) + 4


def interior_cell_count(shape):
    return cell_count(shape) - boundary_cell_count(shape)


def face_count(shape):
    return (shape[0] - 1) * (shape[1] - 2) + (shape[0] - 2) * (shape[1] - 1)


def top_index_iter(shape):
    return xrange(shape[1] * (shape[0] - 1), shape[0] * shape[1])


def bottom_index_iter(shape):
    return xrange(0, shape[1])


def left_index_iter(shape):
    return xrange(0, shape[0] * shape[1], shape[1])


def right_index_iter(shape):
    return xrange(shape[1] - 1, shape[0] * shape[1], shape[1])


def node_xyz(shape, *args):
    """
    Get x, y, and z coordinates for nodes in a structured grid with
    dimensions, *shape*.
    """
    assert(len(shape) == 2)

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

    row_y = np.arange(origin[0], shape[0] * spacing[0] + origin[0],
                         spacing[0])
    col_x = np.arange(origin[1], shape[1] * spacing[1] + origin[0],
                         spacing[1])

    (node_x, node_y) = np.meshgrid(col_x, row_y)
    node_z = np.zeros(node_count)

    node_x.shape = (node_count, )
    node_y.shape = (node_count, )

    return (node_x, node_y, node_z)


def active_cells(shape):
    return np.arange(active_cell_count(shape))


def active_cell_node(shape):
    return cell_node_index(shape)


def node_active_cell(shape):
    node_count = node_count(shape)

    node_ids = np.arange(node_count)
    node_ids.shape = shape

    node_ids[0, :] = _BAD_INDEX_VALUE
    node_ids[:, 0] = _BAD_INDEX_VALUE
    node_ids[-1, :] = _BAD_INDEX_VALUE
    node_ids[:, -1] = _BAD_INDEX_VALUE

    node_ids.shape = (node_count, )

    return node_ids


def cell_node_index(shape):
    node_ids = np.arange(node_count(shape))
    node_ids.shape = shape

    cell_node = node_ids[1:-1, 1:-1].copy()
    cell_node.shape = ((shape[0] - 2) * (shape[1] - 2), )

    return cell_node


def link_lists(shape):
    """
    Link lists:
    For all links, we encode the "from" and "to" nodes, and the face
    (if any) associated with the link. If the link does not intersect a
    face, then face is assigned None.
    For active links, we store the corresponding link ID.

    The numbering scheme for links in RasterModelGrid is illustrated with
    the example of a five-column by four-row grid (each * is a node,
    the lines show links, and the ^ and > symbols indicate the direction
    of each link: up for vertical links, and right for horizontal ones)::

    *--27-->*--28-->*--29-->*--30-->*
    ^       ^       ^       ^       ^
    10      11      12      13      14
    |       |       |       |       |
    *--23-->*--24-->*--25-->*--26-->*
    ^       ^       ^       ^       ^
    5       6       7       8       9   
    |       |       |       |       |
    *--19-->*--20-->*--21-->*--22-->*
    ^       ^       ^       ^       ^
    0       1       2       3       4
    |       |       |       |       |
    *--15-->*--16-->*--17-->*--18-->*
    
    create the fromnode and tonode lists
    """
    (num_rows, num_cols) = shape

    self.link_fromnode = []
    self.link_tonode = []
        
    #   vertical links
    for r in range(0, num_rows-1):
        for c in range(0, num_cols):
            self.link_fromnode.append(c+r*num_cols)
            self.link_tonode.append(c+(r+1)*num_cols)
        
    #   horizontal links
    for r in range(0, num_rows):
        for c in range(0, num_cols-1):
            self.link_fromnode.append(c+r*num_cols)
            self.link_tonode.append(c+r*num_cols+1)
        
    #   convert to np arrays
    self.link_fromnode = np.array(self.link_fromnode)
    self.link_tonode = np.array(self.link_tonode)
    

def node_link_index(shape):
    node_ids = np.arange(np.prod(shape))
    node_ids.shape = shape

    return (from_node_links(node_ids), to_node_links(node_ids))


def to_node_links(node_ids):
    vertical_links = node_ids[1:, :]
    horizontal_links = node_ids[:, 1:]
    return np.concatenate((vertical_links.flat, horizontal_links.flat))


def from_node_links(node_ids):
    vertical_links = node_ids[:-1, :]
    horizontal_links = node_ids[:, :-1]
    return np.concatenate((vertical_links.flat, horizontal_links.flat))


def link_faces(shape, active_links):
    num_links = link_count(shape)

    link_faces = np.empty(num_links, dtype=np.int)

    #active_links = active_links(shape)

    link_faces.fill(_BAD_INDEX_VALUE)
    link_faces[active_links] = np.arange(len(active_links))

    return link_faces


if __name__ == '__main__':
    import doctest
    doctest.testmod()

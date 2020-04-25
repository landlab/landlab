import numpy as np
cimport numpy as np
cimport cython
from libc.stdlib cimport malloc, free


DTYPE = np.int
ctypedef np.int_t DTYPE_t


@cython.boundscheck(False)
def fill_xy_of_node_hex_horizontal(
    shape,
    np.ndarray[np.double_t, ndim=1] x_of_node,
    np.ndarray[np.double_t, ndim=1] y_of_node,
):
    """Get x and y coordinates for each node."""
    cdef long n_rows = shape[0]
    cdef long n_cols = shape[1]
    cdef long n_nodes = len(x_of_node)
    cdef long longest_row = n_rows // 2
    cdef long row
    cdef long offset
    cdef long *size_of_row
    cdef long *offset_to_row
    cdef double x0

    try:
        size_of_row = <long*>malloc(sizeof(long) * n_rows)
        offset_to_row = <long*>malloc(sizeof(long) * (n_rows + 1))

        size_of_row[0] = shape[1]
        for row in range(1, longest_row + 1):
            size_of_row[row] = size_of_row[row - 1] + 1
        for row in range(longest_row + 1, n_rows):
            size_of_row[row] = size_of_row[row - 1] - 1

        offset_to_row[0] = 0
        for row in range(1, n_rows + 1):
            offset_to_row[row] = offset_to_row[row - 1] + size_of_row[row - 1]

        x0 = longest_row / 2.
        for row in range(longest_row + 1):
            for offset in range(offset_to_row[row], offset_to_row[row + 1]):
                x_of_node[offset] = x0 + offset - offset_to_row[row]
                y_of_node[offset] = row
            x0 -= .5

        x0 += 1.
        for row in range(longest_row + 1, n_rows):
            for offset in range(offset_to_row[row], offset_to_row[row + 1]):
                x_of_node[offset] = x0 + offset - offset_to_row[row]
                y_of_node[offset] = row
            x0 += .5
    finally:
        free(offset_to_row)
        free(size_of_row)


@cython.boundscheck(False)
def fill_xy_of_node_hex_vertical(
    shape,
    np.ndarray[np.double_t, ndim=1] x_of_node,
    np.ndarray[np.double_t, ndim=1] y_of_node,
):
    cdef long n_nodes = len(x_of_node)
    cdef int n_cols = shape[1]
    cdef int longest_column = n_cols // 2
    cdef int size_of_longest_column = longest_column + shape[0]
    cdef int n_rows = 2 * size_of_longest_column - 1
    cdef long size_of_odd_row = n_cols // 2
    cdef long size_of_even_row = n_cols - size_of_odd_row
    cdef long n_middle_rows = 2 * shape[0] - 1
    cdef long n_top_rows = n_cols // 2
    cdef long row
    cdef long n
    cdef long *size_of_row
    cdef long *offset_to_row
    cdef float x0
    cdef long offset

    try:
        size_of_row = <long*>malloc(sizeof(long) * n_rows)
        offset_to_row = <long*>malloc(sizeof(long) * (n_rows + 1))

        # Bottom and top rows
        for row in range(n_top_rows):
            size_of_row[row] = row + 1
            size_of_row[n_rows - row - 1] = row + 1

        # Middle rows
        for row in range(0, n_middle_rows, 2):
            size_of_row[row + n_top_rows] = size_of_even_row
            size_of_row[row + 1 + n_top_rows] = size_of_odd_row

        offset_to_row[0] = 0
        for row in range(1, n_rows + 1):
            offset_to_row[row] = offset_to_row[row - 1] + size_of_row[row - 1]

        for row in range(0, n_rows, 2):
            # offset = offset_to_row[row]
            x0 = - (size_of_row[row] - 1) // 2
            # for n in range(size_of_row[row]):
            #     x_of_node[offset + n] = x0 + n
            #     y_of_node[offset + n] = row
            for offset in range(offset_to_row[row], offset_to_row[row + 1]):
                x_of_node[offset] = x0 + offset - offset_to_row[row]
                y_of_node[offset] = row

        for row in range(1, n_rows, 2):
            # offset = offset_to_row[row]
            x0 = - size_of_row[row] // 2 + .5
            # for n in range(size_of_row[row]):
            #     x_of_node[offset + n] = x0 + n
            #     y_of_node[offset + n] = row
            for offset in range(offset_to_row[row], offset_to_row[row + 1]):
                x_of_node[offset] = x0 + offset - offset_to_row[row]
                y_of_node[offset] = row

        x0 = x_of_node[offset_to_row[n_top_rows]]
        for n in range(n_nodes):
            x_of_node[n] -= x0
    finally:
        free(offset_to_row)
        free(size_of_row)


@cython.boundscheck(False)
def fill_xy_of_node_rect_vertical(
    shape,
    np.ndarray[np.double_t, ndim=1] x_of_node,
    np.ndarray[np.double_t, ndim=1] y_of_node,
):
    """Get x and y coordinates for each node."""
    cdef long n_cols = shape[1]
    cdef long size_of_odd_row = n_cols // 2
    cdef long size_of_even_row = n_cols - size_of_odd_row
    cdef long n_rows = 2 * shape[0]
    cdef long row
    cdef long offset
    cdef long n

    # even rows
    offset = 0
    for row in range(0, n_rows, 2):
        for n in range(size_of_even_row):
            x_of_node[offset + n] = n * 2
            y_of_node[offset + n] = row / 2.
        offset += n_cols

    # odd rows
    offset = size_of_even_row
    for row in range(1, n_rows, 2):
        for n in range(size_of_odd_row):
            x_of_node[offset + n] = 2 * n + 1
            y_of_node[offset + n] = row / 2.
        offset += n_cols


@cython.boundscheck(False)
def fill_xy_of_node_rect_horizontal(
    shape,
    np.ndarray[np.double_t, ndim=1] x_of_node,
    np.ndarray[np.double_t, ndim=1] y_of_node,
):
    """Get x and y coordinates for each node."""
    cdef int n_nodes = x_of_node.shape[0]
    cdef int stride = 2 * shape[1]
    cdef int n_cols = shape[1]
    cdef int row
    cdef int offset
    cdef int n

    for offset in range(0, n_nodes, stride):
        row = (offset // stride) * 2
        for n in range(offset, offset + n_cols):
            x_of_node[n] = n - offset
            y_of_node[n] = row

    for offset in range(shape[1], n_nodes, stride):
        row = (offset // stride) * 2 + 1
        for n in range(offset, offset + n_cols):
            x_of_node[n] = n - offset + .5
            y_of_node[n] = row


@cython.boundscheck(False)
def get_xy_of_node(shape,
                   np.ndarray[np.double_t, ndim=1] x_of_node,
                   np.ndarray[np.double_t, ndim=1] y_of_node):
    """Get x and y coordinates for each node."""
    cdef int n_nodes = x_of_node.shape[0]
    cdef int stride = 2 * shape[1] + 1
    cdef int n_cols = shape[1]
    cdef int row
    cdef int offset
    cdef int n

    for offset in range(0, n_nodes, stride):
        row = (offset // stride) * 2
        for n in range(offset, offset + n_cols):
            x_of_node[n] = n - offset + .5
            y_of_node[n] = row

    for offset in range(shape[1], n_nodes, stride):
        row = (offset // stride) * 2 + 1
        for n in range(offset, offset + n_cols + 1):
            x_of_node[n] = n - offset
            y_of_node[n] = row


@cython.boundscheck(False)
def get_nodes_at_link(shape,
                   np.ndarray[DTYPE_t, ndim=2] nodes_at_link):
    """Get nodes at the tail and head of each node."""
    cdef int n_links = nodes_at_link.shape[0]
    cdef int n_short_rows = (shape[0] + 1) // 2
    cdef int n_long_rows = shape[0] // 2
    cdef int n_cols = shape[1]
    cdef int n_rows = shape[0]
    cdef int stride = 2 * n_cols + 1
    cdef int row
    cdef int n
    cdef int link

    link = 0
    for row in range(n_rows):
        node = row * stride // 2
        if row % 2 == 0:
            for n in range(n_cols - 1):
                nodes_at_link[link, 0] = n + node
                nodes_at_link[link, 1] = n + node + 1
                link += 1
        else:
            for n in range(n_cols):
                nodes_at_link[link, 0] = n + node
                nodes_at_link[link, 1] = n + node + 1
                link += 1

        if row == n_rows - 1:
            break

        if row % 2 == 0:
            for n in range(n_cols):
                nodes_at_link[link, 0] = n + node
                nodes_at_link[link, 1] = n + node + n_cols
                link += 1
                nodes_at_link[link, 0] = n + node
                nodes_at_link[link, 1] = n + node + n_cols + 1
                link += 1
        else:
            nodes_at_link[link, 0] = node
            nodes_at_link[link, 1] = node + n_cols + 1
            link += 1
            for n in range(1, n_cols):
                nodes_at_link[link, 0] = n + node
                nodes_at_link[link, 1] = n + node + n_cols
                link += 1
                nodes_at_link[link, 0] = n + node
                nodes_at_link[link, 1] = n + node + n_cols + 1
                link += 1
            nodes_at_link[link, 0] = node + n_cols
            nodes_at_link[link, 1] = node + n_cols + n_cols
            link += 1


@cython.boundscheck(False)
def get_links_at_patch(shape,
                       np.ndarray[DTYPE_t, ndim=2] links_at_patch):
    """Get links that bound each patch."""
    cdef int n_patches = links_at_patch.shape[0]
    cdef int n_short_rows = (shape[0] + 1) // 2
    cdef int n_long_rows = shape[0] // 2
    cdef int n_cols = shape[1]
    cdef int n_rows = shape[0]
    cdef int row
    cdef int n
    cdef int patch
    cdef int link
    cdef int link_at_row

    patch = 0
    link_at_row = n_cols - 1
    for row in range(n_rows - 1):
        if row % 2 == 0:
            link = link_at_row + 2
            n_patches = n_cols - 1
        else:
            link = link_at_row + 1
            n_patches = n_cols

        # upward pointing
        links_at_patch[patch, 0] = link
        links_at_patch[patch, 1] = link - 1
        links_at_patch[patch, 2] = link - (n_cols + 1)
        patch += 1
        for n in range(n_patches - 1):
            links_at_patch[patch, 0] = links_at_patch[patch - 1, 0] + 2
            links_at_patch[patch, 1] = links_at_patch[patch - 1, 1] + 2
            links_at_patch[patch, 2] = links_at_patch[patch - 1, 2] + 1
            patch += 1

        if row % 2 == 0:
            link = link_at_row
            n_patches = n_cols
        else:
            link = link_at_row + 1
            n_patches = n_cols - 1

        # downward pointing
        if row % 2 == 0:
            links_at_patch[patch, 0] = link + 2 * n_cols
        else:
            links_at_patch[patch, 0] = link + 2 * n_cols - 1
        links_at_patch[patch, 1] = link
        links_at_patch[patch, 2] = link + 1
        patch += 1
        for n in range(n_patches - 1):
            links_at_patch[patch, 0] = links_at_patch[patch - 1, 0] + 1
            links_at_patch[patch, 1] = links_at_patch[patch - 1, 1] + 2
            links_at_patch[patch, 2] = links_at_patch[patch - 1, 2] + 2
            patch += 1

        if row % 2 == 0:
            link_at_row += 3 * n_cols
        else:
            link_at_row += 3 * n_cols - 1

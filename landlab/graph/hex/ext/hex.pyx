import numpy as np
cimport numpy as np
cimport cython
from libc.stdlib cimport malloc, free


DTYPE = np.int
ctypedef np.int_t DTYPE_t


@cython.boundscheck(False)
def fill_perimeter_nodes(shape, np.ndarray[DTYPE_t, ndim=1] perimeter_nodes):
    cdef int n_rows = shape[0]
    cdef int n_cols = shape[1]
    cdef int n_nodes = n_rows * n_cols
    cdef int i
    cdef int node

    # Right edge
    i = 0
    for node in range(n_cols - 1, n_nodes - 1, n_cols):
        perimeter_nodes[i] = node
        i += 1

    # Top edge
    for node in range(n_nodes - 1, n_nodes - n_cols, - 1):
        perimeter_nodes[i] = node
        i += 1

    # Left edge
    for node in range((n_rows - 1) * n_cols, 0, - n_cols):
        perimeter_nodes[i] = node
        i += 1

    # Bottom edge
    for node in range(0, n_cols - 1):
        perimeter_nodes[i] = node
        i += 1


@cython.boundscheck(False)
def fill_hex_perimeter_nodes(shape,
                             np.ndarray[DTYPE_t, ndim=1] perimeter_nodes):
    cdef int n_rows = shape[0]
    cdef int n_cols = shape[1]
    cdef int n_bottom_rows = (n_rows + (n_rows + 1) % 2) // 2 + 1
    cdef int i, i0
    cdef int node
    cdef int row
    cdef int * nodes_per_row = <int *>malloc(n_rows * sizeof(int))

    try:
        nodes_per_row[0] = n_cols
        for row in range(1, n_bottom_rows):
            nodes_per_row[row] = nodes_per_row[row - 1] + 1
        for row in range(n_bottom_rows, n_rows):
            nodes_per_row[row] = nodes_per_row[row - 1] - 1

        # Right edge
        perimeter_nodes[0] = n_cols - 1
        row = 1
        for i in range(1, n_rows):
            perimeter_nodes[i] = perimeter_nodes[i - 1] + nodes_per_row[row]
            row += 1

        # Top edge
        i0 = i + 1
        for i in range(i0, i0 + nodes_per_row[n_rows - 1] - 1):
            perimeter_nodes[i] = perimeter_nodes[i - 1] - 1

        # Left edge
        i0 = i + 1
        row = n_rows - 2
        for i in range(i0, i0 + n_rows - 1):
            perimeter_nodes[i] = perimeter_nodes[i - 1] - nodes_per_row[row]
            row -= 1

        # Bottom edge
        i0 = i + 1
        for i in range(i0, i0 + n_cols - 2):
            perimeter_nodes[i] = perimeter_nodes[i - 1] + 1
    finally:
        free(nodes_per_row)


@cython.boundscheck(False)
def fill_hex_xy_of_node(shape,
                        np.ndarray[np.double_t, ndim=1] x_of_node,
                        np.ndarray[np.double_t, ndim=1] y_of_node):
    """Get x and y coordinates for each node."""
    cdef int n_nodes = x_of_node.size
    cdef int n_cols = shape[1]
    cdef int node
    cdef int row
    cdef int col
    cdef int offset
    cdef int n
    cdef double x0

    node = 0
    x0 = 0.
    for row in range(shape[0] / 2 + 1):
        for col in range(n_cols):
            x_of_node[node] = x0 + col
            y_of_node[node] = row
            node += 1
        x0 -= .5
        n_cols += 1

    x0 += 1.
    n_cols -= 2
    for row in range(shape[0] / 2 + 1, shape[0]):
        for col in range(n_cols):
            x_of_node[node] = x0 + col
            y_of_node[node] = row
            node += 1
        x0 += .5
        n_cols -= 1


@cython.boundscheck(False)
def fill_xy_of_node(shape,
                    np.ndarray[np.double_t, ndim=1] x_of_node,
                    np.ndarray[np.double_t, ndim=1] y_of_node):
    """Get x and y coordinates for each node."""
    cdef int n_nodes = x_of_node.size
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
    cdef int n_nodes = x_of_node.size
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

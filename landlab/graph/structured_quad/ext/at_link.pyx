import numpy as np
cimport numpy as np
cimport cython

DTYPE = np.int
ctypedef np.int_t DTYPE_t


@cython.boundscheck(False)
def fill_horizontal_links(shape, np.ndarray[DTYPE_t, ndim=1] horizontal_links):
    cdef int n_rows = shape[0]
    cdef int n_cols = shape[1]
    cdef int n_links = n_rows * (n_cols - 1) + (n_rows - 1) * n_cols
    cdef int link_stride = 2 * n_cols - 1
    cdef int i, n, link

    i = 0
    for link in range(0, n_links, link_stride):
        for n in range(n_cols - 1):
            horizontal_links[i] = link + n
            i += 1


@cython.boundscheck(False)
def fill_vertical_links(shape, np.ndarray[DTYPE_t, ndim=1] vertical_links):
    cdef int n_rows = shape[0]
    cdef int n_cols = shape[1]
    cdef int link_stride = 2 * n_cols - 1
    cdef int n_links = n_rows * (n_cols - 1) + (n_rows - 1) * n_cols
    cdef int i, n, link

    i = 0
    for link in range(n_cols - 1, n_links, link_stride):
        for n in range(n_cols):
            vertical_links[i] = link + n
            i += 1


@cython.boundscheck(False)
def fill_patches_at_link(shape, np.ndarray[DTYPE_t, ndim=2] patches_at_link):
    cdef int link
    cdef int patch
    cdef int n_rows = shape[0]
    cdef int n_cols = shape[1]
    cdef int patches_per_row = n_cols - 1
    cdef int n_links = (2 * n_cols - 1) * n_rows - n_cols

    # Interior horizontal links
    for row in range(1, n_rows - 1):
        patch = patches_per_row * (row - 1)
        link = (2 * n_cols - 1) * row
        for link in range(link, link + n_cols - 1):
            patches_at_link[link, 0] = patch
            patches_at_link[link, 1] = patch + patches_per_row
            patch += 1

    # Interior vertical links
    for row in range(0, n_rows - 1):
        patch = patches_per_row * row
        link = n_cols + (2 * n_cols - 1) * row
        for link in range(link, link + n_cols - 2):
            patches_at_link[link, 0] = patch + 1
            patches_at_link[link, 1] = patch
            patch += 1

    # Left edge
    patch = 0
    for link in range(n_cols - 1, n_links, 2 * n_cols - 1):
        patches_at_link[link, 0] = patch
        patches_at_link[link, 1] = - 1
        patch += patches_per_row

    # Right edge
    patch = patches_per_row - 1
    for link in range(2 * n_cols - 2, n_links, 2 * n_cols - 1):
        patches_at_link[link, 0] = - 1
        patches_at_link[link, 1] = patch
        patch += patches_per_row

    # Bottom edge
    patch = 0
    for link in range(n_cols - 1):
        patches_at_link[link, 0] = - 1
        patches_at_link[link, 1] = patch
        patch += 1

    # Top edge
    patch = patches_per_row * (n_rows - 2)
    for link in range(n_links - (n_cols - 1), n_links):
        patches_at_link[link, 0] = patch
        patches_at_link[link, 1] = - 1
        patch += 1


@cython.boundscheck(False)
def fill_nodes_at_link(shape, np.ndarray[DTYPE_t, ndim=2] nodes_at_link):
    cdef int row, col
    cdef int link
    cdef int node
    cdef int n_rows = shape[0]
    cdef int n_cols = shape[1]
    cdef int n_links = (2 * n_cols - 1) * n_rows - n_cols
    cdef int links_per_row = 2 * n_cols - 1

    # Horizontal links
    for row in range(n_rows):
        node = row * n_cols
        link = row * links_per_row
        for col in range(n_cols - 1):
            nodes_at_link[link, 0] = node
            nodes_at_link[link, 1] = node + 1
            node += 1
            link += 1

    # Vertical links
    for row in range(0, n_rows - 1):
        node = row * n_cols
        link = row * links_per_row + n_cols - 1
        for col in range(n_cols):
            nodes_at_link[link, 0] = node
            nodes_at_link[link, 1] = node + n_cols
            node += 1
            link += 1

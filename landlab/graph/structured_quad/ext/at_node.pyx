cimport cython
from cython.parallel cimport prange

cimport numpy as np


@cython.boundscheck(False)
@cython.wraparound(False)
def fill_perimeter_nodes(
    shape,
    cython.integral[:] perimeter_nodes,
):
    cdef int n_rows = shape[0]
    cdef int n_cols = shape[1]
    cdef int n_nodes = n_rows * n_cols
    cdef int offset_to_top = n_rows - 1
    cdef int offset_to_left = offset_to_top + n_cols - 1
    cdef int offset_to_bottom = offset_to_left + n_rows - 1
    cdef int row
    cdef int col

    for row in prange(n_rows - 1, nogil=True, schedule="static"):
        perimeter_nodes[row] = (row + 1) * n_cols - 1
        perimeter_nodes[offset_to_left + row] = n_nodes - (row + 1) * n_cols

    for col in prange(n_cols - 1, nogil=True, schedule="static"):
        perimeter_nodes[offset_to_top + col] = n_nodes - 1 - col
        perimeter_nodes[offset_to_bottom + col] = col


@cython.boundscheck(False)
@cython.wraparound(False)
def fill_patches_at_node(
    shape,
    cython.integral [:, :] patches_at_face,
):
    cdef int patch
    cdef int node
    cdef int row
    cdef int col
    cdef int n_rows = shape[0]
    cdef int n_cols = shape[1]
    cdef int patches_per_row = n_cols - 1

    # Bottom row
    for col in prange(1, n_cols - 1, nogil=True, schedule="static"):
        node = col
        patch = node

        patches_at_face[node, 0] = node
        patches_at_face[node, 1] = node - 1
        patches_at_face[node, 2] = - 1
        patches_at_face[node, 3] = - 1


    # Interior nodes
    for row in prange(1, n_rows - 1, nogil=True, schedule="static"):
        patch = patches_per_row * row
        node = row * n_cols + 1
        for col in range(1, n_cols - 1):
            patches_at_face[node, 0] = patch + 1
            patches_at_face[node, 1] = patch
            patches_at_face[node, 2] = patch - patches_per_row
            patches_at_face[node, 3] = patch - patches_per_row + 1

            node = node + 1
            patch = patch + 1

    # Top row
    for col in prange(1, n_cols - 1, nogil=True, schedule="static"):
        node = (n_rows - 1) * n_cols + col
        patch = (n_rows - 2) * (n_cols - 1) + col - 1

        patches_at_face[node, 0] = - 1
        patches_at_face[node, 1] = - 1
        patches_at_face[node, 2] = patch
        patches_at_face[node, 3] = patch + 1

    # Left column
    for row in prange(1, n_rows - 1, nogil=True, schedule="static"):
        node = row * n_cols
        patch = row * (n_cols - 1)

        patches_at_face[node, 0] = patch
        patches_at_face[node, 1] = - 1
        patches_at_face[node, 2] = - 1
        patches_at_face[node, 3] = patch - patches_per_row

    # Right column
    for row in prange(1, n_rows - 1, nogil=True, schedule="static"):
        node = row * n_cols + n_cols - 1
        patch = row * (n_cols - 1) + n_cols - 2

        patches_at_face[node, 0] = - 1
        patches_at_face[node, 1] = patch
        patches_at_face[node, 2] = patch - patches_per_row
        patches_at_face[node, 3] = - 1

    # Corners
    node = 0
    patches_at_face[node, 0] = 0
    patches_at_face[node, 1] = - 1
    patches_at_face[node, 2] = - 1
    patches_at_face[node, 3] = - 1

    node = n_cols - 1
    patches_at_face[node, 0] = - 1
    patches_at_face[node, 1] = patches_per_row - 1
    patches_at_face[node, 2] = - 1
    patches_at_face[node, 3] = - 1

    node = (n_rows - 1) * n_cols
    patches_at_face[node, 0] = - 1
    patches_at_face[node, 1] = - 1
    patches_at_face[node, 2] = - 1
    patches_at_face[node, 3] = patches_per_row * (n_rows - 2)

    node = n_rows * n_cols - 1
    patches_at_face[node, 0] = - 1
    patches_at_face[node, 1] = - 1
    patches_at_face[node, 2] = patches_per_row * (n_rows - 1) - 1
    patches_at_face[node, 3] = - 1


@cython.boundscheck(False)
@cython.wraparound(False)
def fill_links_at_node(
    shape,
    cython.integral [:, :] links_at_node,
):
    cdef int n_rows = shape[0]
    cdef int n_cols = shape[1]
    cdef int n_nodes = n_rows * n_cols
    cdef int links_per_row = 2 * n_cols - 1
    cdef int link
    cdef int node
    cdef int row
    cdef int col

    # Bottom nodes
    for col in prange(1, n_cols - 1, nogil=True, schedule="static"):
        node = col
        link = col - 1

        links_at_node[node, 0] = link + 1
        links_at_node[node, 1] = link + n_cols
        links_at_node[node, 2] = link
        links_at_node[node, 3] = - 1


    # Interior nodes
    for row in prange(1, n_rows - 1, nogil=True, schedule="static"):
        link = row * links_per_row
        node = row * n_cols + 1
        for col in range(1, n_cols - 1):
            links_at_node[node, 0] = link + 1
            links_at_node[node, 1] = link + n_cols
            links_at_node[node, 2] = link
            links_at_node[node, 3] = link - (n_cols - 1)

            node = node + 1
            link = link + 1

    # Top nodes
    for col in prange(1, n_cols - 1, nogil=True, schedule="static"):
        node = (n_rows - 1) * n_cols + col
        link = (n_rows - 1) * links_per_row + col - 1

        links_at_node[node, 0] = link + 1
        links_at_node[node, 1] = - 1
        links_at_node[node, 2] = link
        links_at_node[node, 3] = link - (n_cols - 1)


    # Left nodes
    for row in prange(1, n_rows - 1, nogil=True, schedule="static"):
        node = row * n_cols
        link = row * links_per_row - 1

        links_at_node[node, 0] = link + 1
        links_at_node[node, 1] = link + n_cols
        links_at_node[node, 2] = - 1
        links_at_node[node, 3] = link - (n_cols - 1)


    # Right nodes
    for row in prange(1, n_rows - 1, nogil=True, schedule="static"):
        node = row * n_cols + n_cols - 1
        link = row * links_per_row + n_cols - 2
    
        links_at_node[node, 0] = - 1
        links_at_node[node, 1] = link + n_cols
        links_at_node[node, 2] = link
        links_at_node[node, 3] = link - (n_cols - 1)

    node = 0
    links_at_node[node, 0] = 0
    links_at_node[node, 1] = n_cols - 1
    links_at_node[node, 2] = - 1
    links_at_node[node, 3] = - 1

    node = n_cols - 1
    links_at_node[node, 0] = - 1
    links_at_node[node, 1] = 2 * n_cols - 2
    links_at_node[node, 2] = n_cols - 2
    links_at_node[node, 3] = - 1

    node = n_nodes - n_cols
    links_at_node[node, 0] = (n_rows - 1) * links_per_row
    links_at_node[node, 1] = - 1
    links_at_node[node, 2] = - 1
    links_at_node[node, 3] = (n_rows - 1) * links_per_row - n_cols

    node = n_nodes - 1
    links_at_node[node, 0] = - 1
    links_at_node[node, 1] = - 1
    links_at_node[node, 2] = n_rows * links_per_row - n_cols - 1
    links_at_node[node, 3] = (n_rows - 1) * links_per_row - 1


@cython.boundscheck(False)
@cython.wraparound(False)
def fill_link_dirs_at_node(
    shape,
    np.int8_t [:, :] link_dirs_at_node,
):
    cdef int n_rows = shape[0]
    cdef int n_cols = shape[1]
    cdef int n_nodes = n_rows * n_cols
    cdef int node
    cdef int row
    cdef int col

    # Bottom nodes
    for node in prange(1, n_cols - 1, nogil=True, schedule="static"):
        link_dirs_at_node[node, 0] = - 1
        link_dirs_at_node[node, 1] = - 1
        link_dirs_at_node[node, 2] = 1
        link_dirs_at_node[node, 3] = 0

    # Interior nodes
    for row in prange(1, n_rows - 1, nogil=True, schedule="static"):
        node = row * n_cols + 1
        for col in range(1, n_cols - 1):
            link_dirs_at_node[node, 0] = - 1
            link_dirs_at_node[node, 1] = - 1
            link_dirs_at_node[node, 2] = 1
            link_dirs_at_node[node, 3] = 1

            node = node + 1

    # Top nodes
    for node in prange(n_nodes - n_cols + 1, n_nodes - 1, nogil=True, schedule="static"):
        link_dirs_at_node[node, 0] = -1
        link_dirs_at_node[node, 1] = 0
        link_dirs_at_node[node, 2] = 1
        link_dirs_at_node[node, 3] = 1

    # Right nodes
    for row in prange(1, n_rows - 1, nogil=True, schedule="static"):
        node = row * n_cols
        link_dirs_at_node[node, 0] = - 1
        link_dirs_at_node[node, 1] = - 1
        link_dirs_at_node[node, 2] = 0
        link_dirs_at_node[node, 3] = 1

        node = node + n_cols - 1
        link_dirs_at_node[node, 0] = 0
        link_dirs_at_node[node, 1] = - 1
        link_dirs_at_node[node, 2] = 1
        link_dirs_at_node[node, 3] = 1

    node = 0
    link_dirs_at_node[node, 0] = - 1
    link_dirs_at_node[node, 1] = - 1
    link_dirs_at_node[node, 2] = 0
    link_dirs_at_node[node, 3] = 0

    node = n_cols - 1
    link_dirs_at_node[node, 0] = 0
    link_dirs_at_node[node, 1] = - 1
    link_dirs_at_node[node, 2] = 1
    link_dirs_at_node[node, 3] = 0

    node = n_nodes - n_cols
    link_dirs_at_node[node, 0] = - 1
    link_dirs_at_node[node, 1] = 0
    link_dirs_at_node[node, 2] = 0
    link_dirs_at_node[node, 3] = 1

    node = n_nodes - 1
    link_dirs_at_node[node, 0] = 0
    link_dirs_at_node[node, 1] = 0
    link_dirs_at_node[node, 2] = 1
    link_dirs_at_node[node, 3] = 1

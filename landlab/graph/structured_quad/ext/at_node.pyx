cimport cython
from cython.parallel cimport prange
from libc.stdint cimport int8_t

ctypedef fused id_t:
    cython.integral
    long long


@cython.boundscheck(False)
@cython.wraparound(False)
def fill_perimeter_nodes(
    shape,
    id_t [:] perimeter_nodes,
):
    cdef long n_rows = shape[0]
    cdef long n_cols = shape[1]
    cdef long n_nodes = n_rows * n_cols
    cdef long offset_to_top = n_rows - 1
    cdef long offset_to_left = offset_to_top + n_cols - 1
    cdef long offset_to_bottom = offset_to_left + n_rows - 1
    cdef long row
    cdef long col

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
    id_t [:, :] patches_at_node,
):
    cdef long patch
    cdef long node
    cdef long row
    cdef long col
    cdef long n_rows = shape[0]
    cdef long n_cols = shape[1]
    cdef long patches_per_row = n_cols - 1
    cdef long first_node
    cdef long first_patch

    # Bottom row
    for node in prange(n_cols, nogil=True, schedule="static"):
        patch = node

        patches_at_node[node, 0] = patch
        patches_at_node[node, 1] = patch - 1
        patches_at_node[node, 2] = - 1
        patches_at_node[node, 3] = - 1
    patches_at_node[0, 1] = -1
    patches_at_node[n_cols - 1, 0] = -1

    # Interior nodes
    for row in prange(1, n_rows - 1, nogil=True, schedule="static"):
        first_node = row * n_cols
        patch = row * patches_per_row
        for node in range(first_node, first_node + n_cols):
            patches_at_node[node, 0] = patch
            patches_at_node[node, 1] = patch - 1
            patches_at_node[node, 2] = patch - patches_per_row - 1
            patches_at_node[node, 3] = patch - patches_per_row

            patch = patch + 1
        patches_at_node[first_node, 1] = -1
        patches_at_node[first_node, 2] = -1
        patches_at_node[first_node + n_cols - 1, 0] = -1
        patches_at_node[first_node + n_cols - 1, 3] = -1

    # Top row
    first_node = (n_rows - 1) * n_cols
    first_patch = (n_rows - 2) * patches_per_row
    for col in prange(n_cols, nogil=True, schedule="static"):
        node = first_node + col

        patches_at_node[node, 0] = - 1
        patches_at_node[node, 1] = - 1
        patches_at_node[node, 2] = first_patch - 1 + col
        patches_at_node[node, 3] = first_patch + col
    patches_at_node[first_node, 2] = - 1
    patches_at_node[first_node + n_cols - 1, 3] = - 1


@cython.boundscheck(False)
@cython.wraparound(False)
def fill_links_at_node(
    shape,
    id_t [:, :] links_at_node,
):
    cdef long n_rows = shape[0]
    cdef long n_cols = shape[1]
    cdef long vertical_links_per_row = n_cols
    cdef long horizontal_links_per_row = n_cols - 1
    cdef long links_per_row = horizontal_links_per_row + vertical_links_per_row
    cdef long link
    cdef long node
    cdef long row
    cdef long col
    cdef long first_node
    cdef long first_link

    # Bottom row
    for node in prange(n_cols, nogil=True, schedule="static"):
        link = node
        links_at_node[node, 0] = link
        links_at_node[node, 1] = link + horizontal_links_per_row
        links_at_node[node, 2] = link - 1
        links_at_node[node, 3] = - 1
    links_at_node[0, 2] = - 1
    links_at_node[n_cols - 1, 0] = - 1

    # Middle rows
    for row in prange(1, n_rows - 1, nogil=True, schedule="static"):
        first_node = row * n_cols
        link = row * links_per_row
        for node in range(first_node, first_node + n_cols):
            links_at_node[node, 0] = link
            links_at_node[node, 1] = link + horizontal_links_per_row
            links_at_node[node, 2] = link - 1
            links_at_node[node, 3] = link - vertical_links_per_row

            link = link + 1
        links_at_node[first_node, 2] = -1
        links_at_node[first_node + n_cols - 1, 0] = -1

    # Top row
    first_node = (n_rows - 1) * n_cols
    first_link = (n_rows - 1) * links_per_row
    for col in prange(n_cols, nogil=True, schedule="static"):
        node = first_node + col
        link = first_link + col
        links_at_node[node, 0] = link
        links_at_node[node, 1] = - 1
        links_at_node[node, 2] = link - 1
        links_at_node[node, 3] = link - vertical_links_per_row
    links_at_node[first_node, 2] = -1
    links_at_node[first_node + n_cols - 1, 0] = -1


@cython.boundscheck(False)
@cython.wraparound(False)
def fill_link_dirs_at_node(
    shape,
    int8_t [:, :] link_dirs_at_node,
):
    cdef long n_rows = shape[0]
    cdef long n_cols = shape[1]
    cdef long first_node
    cdef long node
    cdef long row

    # Bottom row
    for node in prange(n_cols, nogil=True, schedule="static"):
        link_dirs_at_node[node, 0] = - 1
        link_dirs_at_node[node, 1] = - 1
        link_dirs_at_node[node, 2] = 1
        link_dirs_at_node[node, 3] = 0
    link_dirs_at_node[0, 2] = 0
    link_dirs_at_node[n_cols - 1, 0] = 0

    # Middle rows
    for row in prange(1, n_rows - 1, nogil=True, schedule="static"):
        first_node = row * n_cols

        for node in range(first_node, first_node + n_cols):
            link_dirs_at_node[node, 0] = - 1
            link_dirs_at_node[node, 1] = - 1
            link_dirs_at_node[node, 2] = 1
            link_dirs_at_node[node, 3] = 1
        link_dirs_at_node[first_node, 2] = 0
        link_dirs_at_node[first_node + n_cols - 1, 0] = 0

    # Top row
    first_node = (n_rows - 1) * n_cols
    for node in prange(first_node, first_node + n_cols, nogil=True, schedule="static"):
        link_dirs_at_node[node, 0] = -1
        link_dirs_at_node[node, 1] = 0
        link_dirs_at_node[node, 2] = 1
        link_dirs_at_node[node, 3] = 1
    link_dirs_at_node[first_node, 2] = 0
    link_dirs_at_node[first_node + n_cols - 1, 0] = 0

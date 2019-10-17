import numpy as np
cimport numpy as np
cimport cython

DTYPE = np.int
ctypedef np.int_t DTYPE_t
INT8TYPE = np.int8
ctypedef np.int8_t INT8TYPE_t


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
def fill_patches_at_node(shape, np.ndarray[DTYPE_t, ndim=2] patches_at_face):
    cdef int patch
    cdef int node
    cdef int row
    cdef int col
    cdef int n_rows = shape[0]
    cdef int n_cols = shape[1]
    cdef int patches_per_row = n_cols - 1

    # Bottom row
    patch = 0
    for node in range(1, n_cols - 1):
        patches_at_face[node, 0] = patch + 1
        patches_at_face[node, 1] = patch
        patches_at_face[node, 2] = - 1
        patches_at_face[node, 3] = - 1

        patch += 1

    # Interior nodes
    for row in range(1, n_rows - 1):
        patch = patches_per_row * row
        node = row * n_cols + 1
        for col in range(1, n_cols - 1):
            patches_at_face[node, 0] = patch + 1
            patches_at_face[node, 1] = patch
            patches_at_face[node, 2] = patch - patches_per_row
            patches_at_face[node, 3] = patch - patches_per_row + 1

            node += 1
            patch += 1

    # Top row
    patch = patches_per_row * (n_rows - 2)
    for node in range(n_cols * (n_rows - 1) + 1, n_rows * n_cols - 1):
        patches_at_face[node, 0] = - 1
        patches_at_face[node, 1] = - 1
        patches_at_face[node, 2] = patch
        patches_at_face[node, 3] = patch + 1

        patch += 1

    # Left column
    patch = patches_per_row
    for node in range(n_cols, (n_rows - 1) * n_cols, n_cols):
        patches_at_face[node, 0] = patch
        patches_at_face[node, 1] = - 1
        patches_at_face[node, 2] = - 1
        patches_at_face[node, 3] = patch - patches_per_row

        patch += patches_per_row

    # Right columns
    patch = patches_per_row - 1
    for node in range(2 * n_cols - 1, (n_rows - 1) * n_cols, n_cols):
        patches_at_face[node, 0] = - 1
        patches_at_face[node, 1] = patch + patches_per_row
        patches_at_face[node, 2] = patch
        patches_at_face[node, 3] = - 1

        patch += patches_per_row

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
def fill_links_at_node(shape, np.ndarray[DTYPE_t, ndim=2] links_at_node):
    cdef int n_rows = shape[0]
    cdef int n_cols = shape[1]
    cdef int n_nodes = n_rows * n_cols
    cdef int links_per_row = 2 * n_cols - 1
    cdef int patches_per_row = n_cols - 1

    # Bottom nodes
    link = 0
    for node in range(1, n_cols - 1):
        links_at_node[node, 0] = link + 1
        links_at_node[node, 1] = link + n_cols
        links_at_node[node, 2] = link
        links_at_node[node, 3] = - 1

        link += 1

    # Interior nodes
    for row in range(1, n_rows - 1):
        link = row * links_per_row
        node = row * n_cols + 1
        for col in range(1, n_cols - 1):
            links_at_node[node, 0] = link + 1
            links_at_node[node, 1] = link + n_cols
            links_at_node[node, 2] = link
            links_at_node[node, 3] = link - (n_cols - 1)

            node += 1
            link += 1

    # Top nodes
    link = (n_rows - 1) * links_per_row
    for node in range(n_nodes - n_cols + 1, n_nodes - 1):
        links_at_node[node, 0] = link + 1
        links_at_node[node, 1] = - 1
        links_at_node[node, 2] = link
        links_at_node[node, 3] = link - (n_cols - 1)

        link += 1

    # Left nodes
    link = links_per_row - 1
    for node in range(n_cols, n_nodes, n_cols):
        links_at_node[node, 0] = link + 1
        links_at_node[node, 1] = link + n_cols
        links_at_node[node, 2] = - 1
        links_at_node[node, 3] = link - (n_cols - 1)

        link += links_per_row

    # Right nodes
    link = links_per_row + (n_cols - 2)
    for node in range(2 * n_cols - 1, n_nodes, n_cols):
        links_at_node[node, 0] = - 1
        links_at_node[node, 1] = link + n_cols
        links_at_node[node, 2] = link
        links_at_node[node, 3] = link - (n_cols - 1)

        link += links_per_row

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
def fill_link_dirs_at_node(shape,
                           np.ndarray[INT8TYPE_t, ndim=2] link_dirs_at_node):
    cdef int n_rows = shape[0]
    cdef int n_cols = shape[1]
    cdef int n_nodes = n_rows * n_cols
    cdef int links_per_row = 2 * n_cols - 1
    cdef int patches_per_row = n_cols - 1

    # Bottom nodes
    for node in range(1, n_cols - 1):
        link_dirs_at_node[node, 0] = - 1
        link_dirs_at_node[node, 1] = - 1
        link_dirs_at_node[node, 2] = 1
        link_dirs_at_node[node, 3] = 0

    # Interior nodes
    for row in range(1, n_rows - 1):
        node = row * n_cols + 1
        for col in range(1, n_cols - 1):
            link_dirs_at_node[node, 0] = - 1
            link_dirs_at_node[node, 1] = - 1
            link_dirs_at_node[node, 2] = 1
            link_dirs_at_node[node, 3] = 1

            node += 1

    # Top nodes
    for node in range(n_nodes - n_cols + 1, n_nodes - 1):
        link_dirs_at_node[node, 0] = -1
        link_dirs_at_node[node, 1] = 0
        link_dirs_at_node[node, 2] = 1
        link_dirs_at_node[node, 3] = 1

    # Left nodes
    for node in range(n_cols, n_nodes, n_cols):
        link_dirs_at_node[node, 0] = - 1
        link_dirs_at_node[node, 1] = - 1
        link_dirs_at_node[node, 2] = 0
        link_dirs_at_node[node, 3] = 1

    # Right nodes
    for node in range(2 * n_cols - 1, n_nodes, n_cols):
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

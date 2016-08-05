import numpy as np
cimport numpy as np
cimport cython

DTYPE = np.int
ctypedef np.int_t DTYPE_t

@cython.boundscheck(False)
def get_node_at_cell(shape, np.ndarray[DTYPE_t, ndim=1] node_at_cell):
    """Get node contained in a cell.

    Parameters
    ----------
    shape : tuple of int
        Shape of the grid as `(n_rows, n_cols)`.
    node_at_cell : ndarray of int
        Buffer into which to place node identifiers.
    """
    cdef int cell
    cdef int cell_rows = shape[0] - 2
    cdef int cell_cols = shape[1] - 2
    cdef int node_cols = shape[1]
    cdef int row_offset
    cdef int row
    cdef int col

    cell = 0
    row_offset = shape[1] + 1
    for row in range(cell_rows):
        for col in range(cell_cols):
            node_at_cell[cell] = row_offset + col
            cell += 1
        row_offset += node_cols


@cython.boundscheck(False)
def get_nodes_at_face(shape, np.ndarray[DTYPE_t, ndim=2] nodes_at_face):
    """Get nodes on either side of a face.

    Parameters
    ----------
    shape : tuple of int
        Shape of the grid as `(n_rows, n_cols)`.
    nodes_at_face : ndarray of int, shape `(n_faces, 2)`
        Buffer into which to place node identifiers.
    """
    cdef int face
    cdef int node
    cdef int row
    cdef int col
    cdef int n_rows = shape[0]
    cdef int n_cols = shape[1]

    # Horizontal faces first
    face = 0
    for row in range(n_rows - 1):
        for col in range(1, n_cols - 1):
            node = row * n_cols + col
            nodes_at_face[face, 0] = node
            nodes_at_face[face, 1] = node + n_cols
            face += 1
        face += n_cols - 1

    # Vertical faces next
    face = n_cols - 2
    for row in range(1, n_rows - 1):
        for col in range(n_cols - 1):
            node = row * n_cols + col
            nodes_at_face[face, 0] = node
            nodes_at_face[face, 1] = node + 1
            face += 1
        face += n_cols - 2


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


@cython.boundscheck(False)
def fill_links_at_patch(shape, np.ndarray[DTYPE_t, ndim=2] links_at_patch):
    cdef int n_rows = shape[0]
    cdef int n_cols = shape[1]
    cdef int links_per_row = 2 * n_cols - 1
    cdef int patches_per_row = n_cols - 1

    for row in range(n_rows - 1):
        link = row * links_per_row + n_cols
        patch = row * patches_per_row
        for col in range(n_cols - 1):
            links_at_patch[patch, 0] = link
            links_at_patch[patch, 1] = link + n_cols - 1
            links_at_patch[patch, 2] = link - 1
            links_at_patch[patch, 3] = link - n_cols

            patch += 1
            link += 1


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
                           np.ndarray[DTYPE_t, ndim=2] link_dirs_at_node):
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

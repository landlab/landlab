cimport cython
from cython.parallel cimport prange


@cython.boundscheck(False)
@cython.wraparound(False)
def fill_horizontal_links(
    shape,
    cython.integral [:] horizontal_links,
):
    cdef int n_rows = shape[0]
    cdef int n_cols = shape[1]
    cdef int link_stride = 2 * n_cols - 1
    cdef int i
    cdef int n
    cdef int link
    cdef int row

    for row in prange(n_rows, nogil=True, schedule="static"):
        link = row * link_stride
        i = row * (n_cols - 1)
        for n in range(n_cols - 1):
            horizontal_links[i + n] = link + n


@cython.wraparound(False)
@cython.boundscheck(False)
def fill_vertical_links(
    shape,
    cython.integral[:] vertical_links,
):
    cdef int n_rows = shape[0]
    cdef int n_cols = shape[1]
    cdef int link_stride = 2 * n_cols - 1
    # cdef int n_links = n_rows * (n_cols - 1) + (n_rows - 1) * n_cols
    cdef int i
    cdef int n
    cdef int link
    cdef int row

    for row in prange(n_rows - 1, nogil=True, schedule="static"):
        link = row * link_stride + n_cols - 1
        i = row * n_cols
        for n in range(n_cols):
            vertical_links[i + n] = link + n


@cython.wraparound(False)
@cython.boundscheck(False)
def fill_patches_at_link(
    shape,
    cython.integral[:, :] patches_at_link,
):
    cdef int link
    cdef int patch
    cdef int n_rows = shape[0]
    cdef int n_cols = shape[1]
    cdef int patches_per_row = n_cols - 1
    cdef int link_stride = 2 * n_cols - 1
    cdef int n_links = (2 * n_cols - 1) * n_rows - n_cols
    cdef int row

    # Interior horizontal links
    for row in prange(1, n_rows - 1, nogil=True, schedule="static"):
        patch = patches_per_row * (row - 1)
        link = (2 * n_cols - 1) * row
        for link in range(link, link + n_cols - 1):
            patches_at_link[link, 0] = patch
            patches_at_link[link, 1] = patch + patches_per_row
            patch = patch + 1

    # Interior vertical links
    for row in prange(n_rows - 1, nogil=True, schedule="static"):
        patch = patches_per_row * row
        # link = n_cols + (2 * n_cols - 1) * row
        link = row * link_stride + n_cols
        for link in range(link, link + n_cols - 2):
            patches_at_link[link, 0] = patch + 1
            patches_at_link[link, 1] = patch
            patch = patch + 1

    # Left and right edges
    for row in prange(n_rows - 1, nogil=True, schedule="static"):
        patch = row * patches_per_row

        link = row * link_stride + n_cols - 1
        patches_at_link[link, 0] = patch
        patches_at_link[link, 1] = - 1

        link = link + n_cols - 1
        patches_at_link[link, 0] = - 1
        patches_at_link[link, 1] = patch + n_cols - 2

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


@cython.wraparound(False)
@cython.boundscheck(False)
def fill_nodes_at_link(
    shape,
    cython.integral[:, :] nodes_at_link,
):
    cdef int row, col
    cdef int link
    cdef int node
    cdef int n_rows = shape[0]
    cdef int n_cols = shape[1]
    cdef int links_per_row = 2 * n_cols - 1

    # Horizontal links
    for row in prange(n_rows, nogil=True, schedule="static"):
        node = row * n_cols
        link = row * links_per_row
        for col in range(n_cols - 1):
            nodes_at_link[link, 0] = node
            nodes_at_link[link, 1] = node + 1
            node = node + 1
            link = link + 1

    # Vertical links
    for row in prange(n_rows - 1, nogil=True, schedule="static"):
        node = row * n_cols
        link = row * links_per_row + n_cols - 1
        for col in range(n_cols):
            nodes_at_link[link, 0] = node
            nodes_at_link[link, 1] = node + n_cols
            node = node + 1
            link = link + 1

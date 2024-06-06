cimport cython
from cython.parallel cimport prange

ctypedef fused id_t:
    cython.integral
    long long


@cython.boundscheck(False)
@cython.wraparound(False)
def fill_horizontal_links(
    shape,
    id_t [:] horizontal_links,
):
    cdef long n_rows = shape[0]
    cdef long n_cols = shape[1]
    cdef long horizontal_links_per_row = n_cols - 1
    cdef long vertical_links_per_row = n_cols
    cdef long links_per_row = horizontal_links_per_row + vertical_links_per_row
    cdef long horizontal_link
    cdef long link
    cdef long row
    cdef long col

    for row in prange(n_rows, nogil=True, schedule="static"):
        link = row * links_per_row
        horizontal_link = row * horizontal_links_per_row

        for col in range(horizontal_links_per_row):
            horizontal_links[horizontal_link + col] = link + col


@cython.wraparound(False)
@cython.boundscheck(False)
def fill_vertical_links(
    shape,
    id_t [:] vertical_links,
):
    cdef long n_rows = shape[0]
    cdef long n_cols = shape[1]
    cdef long horizontal_links_per_row = n_cols - 1
    cdef long vertical_links_per_row = n_cols
    cdef long links_per_row = horizontal_links_per_row + vertical_links_per_row
    cdef long link
    cdef long vertical_link
    cdef long row
    cdef long col

    for row in prange(n_rows - 1, nogil=True, schedule="static"):
        link = row * links_per_row + horizontal_links_per_row
        vertical_link = row * vertical_links_per_row

        for col in range(n_cols):
            vertical_links[vertical_link + col] = link + col


@cython.wraparound(False)
@cython.boundscheck(False)
def fill_patches_at_link(
    shape,
    id_t [:, :] patches_at_link,
):
    cdef long n_rows = shape[0]
    cdef long n_cols = shape[1]
    cdef long patches_per_row = n_cols - 1
    cdef long horizontal_links_per_row = n_cols - 1
    cdef long vertical_links_per_row = n_cols
    cdef long links_per_row = horizontal_links_per_row + vertical_links_per_row
    cdef long link
    cdef long patch
    cdef long row
    cdef long col
    cdef long first_link

    for row in prange(n_rows - 1, nogil=True, schedule="static"):
        first_link = row * links_per_row
        patch = row * patches_per_row
        for link in range(first_link, first_link + horizontal_links_per_row):
            patches_at_link[link, 0] = patch - patches_per_row
            patches_at_link[link, 1] = patch
            patch = patch + 1

        first_link = row * links_per_row + horizontal_links_per_row
        patch = row * patches_per_row
        for link in range(first_link, first_link + vertical_links_per_row):
            patches_at_link[link, 0] = patch
            patches_at_link[link, 1] = patch - 1
            patch = patch + 1
        patches_at_link[first_link, 1] = - 1
        patches_at_link[first_link + vertical_links_per_row - 1, 0] = - 1

    # Bottom row
    for link in prange(horizontal_links_per_row, nogil=True, schedule="static"):
        patches_at_link[link, 0] = -1

    # Top row
    first_link = (n_rows - 1) * links_per_row
    patch = (n_rows - 2) * patches_per_row
    for col in prange(horizontal_links_per_row, nogil=True, schedule="static"):
        patches_at_link[first_link + col, 0] = patch + col
        patches_at_link[first_link + col, 1] = -1


@cython.wraparound(False)
@cython.boundscheck(False)
def fill_nodes_at_link(
    shape,
    id_t [:, :] nodes_at_link,
):
    cdef long row
    cdef long col
    cdef long link
    cdef long node
    cdef long n_rows = shape[0]
    cdef long n_cols = shape[1]
    cdef long horizontal_links_per_row = n_cols - 1
    cdef long vertical_links_per_row = n_cols
    cdef long links_per_row = horizontal_links_per_row + vertical_links_per_row

    for row in prange(n_rows - 1, nogil=True, schedule="static"):
        node = row * n_cols
        link = row * links_per_row

        for col in range(horizontal_links_per_row):
            nodes_at_link[link, 0] = node
            nodes_at_link[link, 1] = node + 1
            node = node + 1
            link = link + 1

        node = row * n_cols
        for col in range(vertical_links_per_row):
            nodes_at_link[link, 0] = node
            nodes_at_link[link, 1] = node + n_cols
            node = node + 1
            link = link + 1

    node = (n_rows - 1) * n_cols
    link = (n_rows - 1) * links_per_row
    for col in prange(horizontal_links_per_row, nogil=True, schedule="static"):
        nodes_at_link[link + col, 0] = node + col
        nodes_at_link[link + col, 1] = node + col + 1

cimport cython
from libc.stdlib cimport free
from libc.stdlib cimport malloc

ctypedef fused id_t:
    cython.integral
    long long


@cython.boundscheck(False)
@cython.wraparound(False)
def fill_perimeter_nodes_rect_horizontal(
    shape,
    id_t [:] perimeter_nodes,
):
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
@cython.wraparound(False)
def fill_perimeter_nodes_rect_vertical(
    shape,
    id_t [:] perimeter_nodes,
):
    cdef int n_rows = shape[0]
    cdef int n_cols = shape[1]
    cdef int i
    cdef int offset_to_right_edge
    cdef int offset_to_top_edge
    cdef int offset_to_left_edge
    cdef int offset_to_bottom_edge

    offset_to_top_edge = n_rows - 1
    offset_to_left_edge = offset_to_top_edge + n_cols - 1
    offset_to_bottom_edge = offset_to_left_edge + n_rows - 1
    offset_to_right_edge = offset_to_bottom_edge + n_cols - 1

    if n_cols % 2 == 0:
        perimeter_nodes[0] = n_cols - 1
    else:
        perimeter_nodes[0] = n_cols // 2

    # Right edge
    for i in range(1, offset_to_top_edge + 1):
        perimeter_nodes[i] = perimeter_nodes[i - 1] + n_cols

    # Top edge
    perimeter_nodes[offset_to_top_edge + 1] = perimeter_nodes[offset_to_top_edge]
    if n_cols % 2 == 0:
        perimeter_nodes[offset_to_top_edge + 1] -= n_cols // 2
    else:
        perimeter_nodes[offset_to_top_edge + 1] += n_cols // 2
    for i in range(offset_to_top_edge + 2, offset_to_left_edge + 1, 2):
        perimeter_nodes[i] = perimeter_nodes[i - 2] - 1
        perimeter_nodes[i + 1] = perimeter_nodes[i - 1] - 1

    # Left edge
    for i in range(offset_to_left_edge + 1, offset_to_bottom_edge + 1):
        perimeter_nodes[i] = perimeter_nodes[i - 1] - n_cols

    # Bottom edge
    perimeter_nodes[offset_to_bottom_edge] = 0
    for i in range(offset_to_bottom_edge + 2, offset_to_right_edge, 2):
        perimeter_nodes[i] = perimeter_nodes[i - 2] + 1

    perimeter_nodes[offset_to_bottom_edge + 1] = (n_cols + 1) // 2
    for i in range(offset_to_bottom_edge + 3, offset_to_right_edge, 2):
        perimeter_nodes[i] = perimeter_nodes[i - 2] + 1


@cython.boundscheck(False)
@cython.wraparound(False)
def fill_perimeter_nodes_hex_horizontal(
    shape,
    id_t [:] perimeter_nodes,
):
    cdef int n_rows = shape[0]
    cdef int n_cols = shape[1]
    cdef int longest_row = n_rows // 2
    cdef int offset_to_right_edge
    cdef int offset_to_top_edge
    cdef int offset_to_left_edge
    cdef int offset_to_bottom_edge
    cdef int i
    cdef int row
    cdef int * nodes_per_row

    try:
        nodes_per_row = <int *>malloc(n_rows * sizeof(int))

        nodes_per_row[0] = n_cols
        for row in range(1, longest_row + 1):
            nodes_per_row[row] = nodes_per_row[row - 1] + 1
        for row in range(longest_row + 1, n_rows):
            nodes_per_row[row] = nodes_per_row[row - 1] - 1

        offset_to_top_edge = n_rows - 1
        offset_to_left_edge = offset_to_top_edge + nodes_per_row[n_rows - 1] - 1
        offset_to_bottom_edge = offset_to_left_edge + n_rows - 1
        offset_to_right_edge = offset_to_bottom_edge + n_cols - 1

        perimeter_nodes[0] = n_cols - 1

        # Right edge
        row = 1
        for i in range(1, offset_to_top_edge + 1):
            perimeter_nodes[i] = perimeter_nodes[i - 1] + nodes_per_row[row]
            row += 1

        # Top edge
        for i in range(offset_to_top_edge + 1, offset_to_left_edge + 1):
            perimeter_nodes[i] = perimeter_nodes[i - 1] - 1

        # Left edge
        row = n_rows - 2
        for i in range(offset_to_left_edge + 1, offset_to_bottom_edge + 1):
            perimeter_nodes[i] = perimeter_nodes[i - 1] - nodes_per_row[row]
            row -= 1

        # Bottom edge
        for i in range(offset_to_bottom_edge + 1, offset_to_right_edge):
            perimeter_nodes[i] = perimeter_nodes[i - 1] + 1
    finally:
        free(nodes_per_row)


@cython.boundscheck(False)
@cython.wraparound(False)
def fill_perimeter_nodes_hex_vertical(
    shape,
    id_t [:] perimeter_nodes,
):
    cdef int n_rows = shape[0]
    cdef int n_cols = shape[1]
    cdef int longest_col = n_cols // 2
    cdef int max_nodes_per_row = n_cols - n_cols // 2
    cdef int offset_to_right_edge
    cdef int offset_to_top_edge
    cdef int offset_to_left_edge
    cdef int offset_to_bottom_edge
    cdef int i

    if n_cols % 2 == 1:
        offset_to_top_edge = n_rows - 1
    else:
        offset_to_top_edge = n_rows
    offset_to_left_edge = offset_to_top_edge + n_cols - 1
    offset_to_bottom_edge = offset_to_left_edge + n_rows - 1
    offset_to_right_edge = offset_to_bottom_edge + n_cols - 1

    # Bottom edge
    nodes_per_row = 1
    perimeter_nodes[offset_to_bottom_edge + longest_col] = 0
    for i in range(offset_to_bottom_edge + longest_col + 1, offset_to_right_edge):
        nodes_per_row += 1
        perimeter_nodes[i] = perimeter_nodes[i - 1] + nodes_per_row

    perimeter_nodes[0] = perimeter_nodes[offset_to_right_edge - 1] + max_nodes_per_row

    # Right edge
    for i in range(1, offset_to_top_edge + 1):
        perimeter_nodes[i] = perimeter_nodes[i - 1] + n_cols

    # Top edge
    nodes_per_row = max_nodes_per_row
    for i in range(offset_to_top_edge + 1, offset_to_left_edge - longest_col + 1):
        nodes_per_row -= 1
        perimeter_nodes[i] = perimeter_nodes[i - 1] + nodes_per_row
    nodes_per_row = 1
    for i in range(offset_to_left_edge - longest_col + 1, offset_to_left_edge):
        nodes_per_row += 1
        perimeter_nodes[i] = perimeter_nodes[i - 1] - nodes_per_row
    perimeter_nodes[offset_to_left_edge] = (
        perimeter_nodes[offset_to_left_edge - 1] - max_nodes_per_row
    )

    # Left edge
    for i in range(offset_to_left_edge + 1, offset_to_bottom_edge + 1):
        perimeter_nodes[i] = perimeter_nodes[i - 1] - n_cols

    # Bottom edge
    perimeter_nodes[offset_to_bottom_edge + 1] = (
        perimeter_nodes[offset_to_bottom_edge] - nodes_per_row
    )
    nodes_per_row = max_nodes_per_row
    for i in range(offset_to_bottom_edge + 2, offset_to_bottom_edge + longest_col):
        nodes_per_row -= 1
        perimeter_nodes[i] = perimeter_nodes[i - 1] - nodes_per_row

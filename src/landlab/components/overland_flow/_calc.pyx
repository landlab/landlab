import numpy as np

cimport cython
cimport numpy as np
cimport openmp as omp
from cython.parallel cimport parallel
from cython.parallel cimport prange
from libc.math cimport fabs
from libc.math cimport powf
from libc.math cimport sqrt
from libc.stdint cimport int8_t

ctypedef fused id_t:
    cython.integral
    long long


@cython.boundscheck(False)
def neighbors_at_link(
    id_t [:] links,
    shape,
    id_t [:, :] out,
):
    cdef int stride
    cdef int n_links
    cdef int link
    cdef int i
    cdef bint is_top, is_bottom, is_left, is_right

    stride = 2 * shape[1] - 1
    n_links = (shape[0] - 1) * shape[1] + shape[0] * (shape[1] - 1)

    for i in range(links.shape[0]):
        link = links[i]

        is_top = link > (n_links - stride)
        is_bottom = link < stride
        is_left = link % stride == 0 or (link + shape[1]) % stride == 0
        is_right = (link - (shape[1] - 2)) % stride == 0 or (link + 1) % stride == 0

        if not is_right:
            out[i, 0] = link + 1

        if not is_top:
            out[i, 1] = link + stride

        if not is_left:
            out[i, 2] = link - 1

        if not is_bottom:
            out[i, 3] = link - stride


@cython.boundscheck(False)
@cython.wraparound(False)
def find_max_water_depth(
    cython.floating [:] h_at_node,
    const id_t [:] nodes,
):
    cdef long n_nodes = len(nodes)
    cdef double max_val = 0.0
    cdef double max_local
    cdef double *ptr_max_val
    cdef double *ptr_max_local
    cdef long i
    cdef long node
    cdef omp.omp_lock_t mutex

    ptr_max_val = &max_val
    omp.omp_init_lock(&mutex)

    with nogil, parallel():
        max_local = 0.0
        ptr_max_local = &max_local

        for i in prange(n_nodes, schedule="static"):
            node = nodes[i]
            if h_at_node[node] > ptr_max_local[0]:
                ptr_max_local[0] = h_at_node[node]

        omp.omp_set_lock(&mutex)
        ptr_max_val[0] = max(ptr_max_val[0], max_local)
        omp.omp_unset_lock(&mutex)

    omp.omp_destroy_lock(&mutex)

    return max_val


@cython.boundscheck(False)
@cython.wraparound(False)
@cython.cdivision(True)
def calc_grad_at_some_links(
    const cython.floating [:] value_at_node,
    const id_t [:, :] nodes_at_link,
    const cython.floating [:] length_of_link,
    const id_t [:] links,
    cython.floating [:] out_at_link,
):
    cdef long n_links = len(links)
    cdef long i
    cdef long link
    cdef long head
    cdef long tail

    for i in prange(n_links, nogil=True, schedule="static"):
        link = links[i]
        tail = nodes_at_link[link, 0]
        head = nodes_at_link[link, 1]

        out_at_link[link] = (
            value_at_node[head] - value_at_node[tail]
        ) / length_of_link[link]


@cython.boundscheck(False)
@cython.wraparound(False)
@cython.cdivision(True)
def calc_discharge_at_link(
    shape,
    cython.floating [:] q_at_link,
    const cython.floating [:] h_at_link,
    const cython.floating [:] water_slope_at_link,
    const cython.floating [:] mannings_at_link,
    const double theta,
    const double g,
    const double dt,
):
    cdef long n_rows = shape[0]
    cdef long n_cols = shape[1]
    cdef long n_links = (n_cols - 1) * n_rows + n_cols * (n_rows - 1)
    cdef long link
    cdef np.ndarray[cython.floating, ndim=1] q_mean_at_link = np.ones_like(q_at_link)

    weighted_mean_of_parallel_links(
        shape,
        theta,
        q_at_link,
        q_mean_at_link,
    )

    for link in prange(n_links, nogil=True, schedule="static"):
        q_at_link[link] = powf(h_at_link[link], 7.0 / 3.0) * (
            q_mean_at_link[link]
            - g * dt * h_at_link[link] * water_slope_at_link[link]
        ) / (
            powf(h_at_link[link], 7.0 / 3.0)
            + g * dt * mannings_at_link[link] ** 2 * fabs(q_at_link[link])
        )


@cython.boundscheck(False)
@cython.wraparound(False)
@cython.cdivision(True)
def calc_discharge_at_some_links(
    cython.floating [:] q_at_link,
    const cython.floating [:] q_mean_at_link,
    const cython.floating [:] h_at_link,
    const cython.floating [:] water_slope_at_link,
    const cython.floating [:] mannings_at_link,
    const id_t [:] links,
    const double g,
    const double dt,
):
    cdef long n_links = len(links)
    cdef long link
    cdef long i
    cdef double h_to_seven_thirds
    cdef double numerator
    cdef double denominator
    cdef double g_times_dt = g * dt
    cdef double seven_thirds = 7.0 / 3.0
    cdef double q_manning
    cdef double t
    cdef double w

    for i in prange(n_links, nogil=True, schedule="static"):
        link = links[i]

        if h_at_link[link] > 0.0:
            h_to_seven_thirds = powf(h_at_link[link], seven_thirds)

            t = h_to_seven_thirds / (g * mannings_at_link[link]**2 * q_at_link[link])
            w = t / (t + dt)

            q_manning =  (
                1.0
                / mannings_at_link[link]
                * powf(h_at_link[link], 5.0 / 3.0)
                * powf(water_slope_at_link[link], 0.5)
            )

            numerator = (
                q_mean_at_link[link]
                - g_times_dt * h_at_link[link] * water_slope_at_link[link]
            )

            denominator = (
                1.0
                + g_times_dt * mannings_at_link[link] ** 2 * fabs(q_at_link[link])
                / h_to_seven_thirds
            )

            q_at_link[link] = numerator / denominator * w + (1.0 - w) * q_manning


@cython.boundscheck(False)
@cython.wraparound(False)
def sum_parallel_links(
    cython.numeric[:] out,
    const cython.numeric[:] value_at_link,
    shape,
):
    cdef int n_rows = shape[0]
    cdef int n_cols = shape[1]
    cdef int links_per_row = 2 * shape[1] - 1
    cdef int row, col
    cdef int link

    for row in prange(n_rows, nogil=True, schedule="static"):
        link = row * links_per_row + 1
        for col in range(1, n_cols - 2):
            out[link] = value_at_link[link - 1] + value_at_link[link + 1]
            link = link + 1

    for row in prange(1, n_rows - 2, nogil=True, schedule="static"):
        link = row * links_per_row + n_cols - 1
        for col in range(n_cols):
            out[link] = (
                value_at_link[link - links_per_row]
                + value_at_link[link + links_per_row]
            )
            link = link + 1


@cython.boundscheck(False)
@cython.wraparound(False)
def weighted_mean_of_parallel_links(
    shape,
    const double weight,
    const cython.floating [:] value_at_link,
    cython.floating [:] out,
):
    cdef long n_rows = shape[0]
    cdef long n_cols = shape[1]
    cdef long horizontal_links_per_row = n_cols - 1
    cdef long vertical_links_per_row = n_cols
    cdef long links_per_row = horizontal_links_per_row + vertical_links_per_row
    cdef long row
    cdef long first_link
    cdef long link

    for row in prange(0, n_rows, nogil=True, schedule="static"):
        first_link = links_per_row * row

        link = first_link
        out[link] = _calc_weighted_mean(
            0.0,
            value_at_link[link],
            value_at_link[link + 1],
            weight,
        )
        for link in range(first_link + 1, first_link + horizontal_links_per_row - 1):
            out[link] = _calc_weighted_mean(
                value_at_link[link - 1],
                value_at_link[link],
                value_at_link[link + 1],
                weight,
            )
        link = first_link + horizontal_links_per_row - 1
        out[link] = _calc_weighted_mean(
            value_at_link[link - 1],
            value_at_link[link],
            0.0,
            weight,
        )

    with nogil:
        first_link = horizontal_links_per_row
        for link in range(first_link, first_link + vertical_links_per_row):
            out[link] = _calc_weighted_mean(
                0.0,
                value_at_link[link],
                value_at_link[link + links_per_row],
                weight,
            )
    for row in prange(1, n_rows - 1, nogil=True, schedule="static"):
        first_link = links_per_row * row + horizontal_links_per_row
        for link in range(first_link, first_link + vertical_links_per_row):
            out[link] = _calc_weighted_mean(
                value_at_link[link - links_per_row],
                value_at_link[link],
                value_at_link[link + links_per_row],
                weight,
            )
    with nogil:
        first_link = links_per_row * (n_rows - 2) + horizontal_links_per_row
        for link in range(first_link, first_link + vertical_links_per_row):
            out[link] = _calc_weighted_mean(
                value_at_link[link - links_per_row],
                value_at_link[link],
                0.0,
                weight,
            )


cdef cython.floating _calc_weighted_mean(
    cython.floating value_at_left,
    cython.floating value_at_center,
    cython.floating value_at_right,
    cython.floating weight,
) noexcept nogil:
    return weight * value_at_center + (1.0 - weight) * 0.5 * (
        value_at_left + value_at_right
    )


@cython.boundscheck(False)
@cython.wraparound(False)
def calc_bates_flow_height_at_some_links(
    const cython.floating [:] z_at_node,
    const cython.floating [:] h_at_node,
    const id_t [:, :] nodes_at_link,
    const id_t [:] links,
    cython.floating [:] out_at_link,
):
    """
    Per Bates et al., 2010, this solution needs to find difference
    between the highest water surface in the two cells and the
    highest bed elevation
    """
    cdef long n_links = len(links)
    cdef long i
    cdef long link
    cdef long head
    cdef long tail

    for i in prange(n_links, nogil=True, schedule="static"):
        link = links[i]
        tail = nodes_at_link[link, 0]
        head = nodes_at_link[link, 1]

        out_at_link[link] = max(
            h_at_node[tail] + z_at_node[tail],
            h_at_node[head] + z_at_node[head],
        ) - max(z_at_node[tail], z_at_node[head])


@cython.boundscheck(False)
@cython.wraparound(False)
def adjust_discharge_for_dry_links(
    const cython.floating [:] h_at_link,
    cython.floating [:] q_at_link,
    const id_t [:] links,
):
    cdef long n_links = len(links)
    cdef long link
    cdef long i

    for i in prange(n_links, nogil=True, schedule="static"):
        link = links[i]
        if h_at_link[link] <= 0.0:
            q_at_link[link] = 0.0


@cython.boundscheck(False)
@cython.wraparound(False)
def adjust_supercritical_discharge(
    cython.floating [:] q_at_link,
    const cython.floating [:] h_at_link,
    const id_t [:] links,
    const double g,
    const double froude,
):
    cdef long n_links = len(links)
    cdef long i
    cdef long link
    cdef double root_g = np.sqrt(g)
    cdef double c

    for i in prange(n_links, nogil=True, schedule="static"):
        link = links[i]

        c = h_at_link[link] * sqrt(h_at_link[link]) * root_g

        if fabs(q_at_link[link]) > froude * c:
            if q_at_link[link] < 0.0:
                q_at_link[link] = -c * froude
            else:
                q_at_link[link] = c * froude


@cython.boundscheck(False)
@cython.wraparound(False)
@cython.cdivision(True)
def adjust_unstable_discharge(
    cython.floating [:] q_at_link,
    const cython.floating [:] h_at_link,
    const id_t [:] links,
    const double dx,
    const double dt,
):
    cdef long n_links = len(links)
    cdef long i
    cdef long link
    cdef double dx_over_dt = dx / dt
    cdef double four_dt_over_dx = 4.0 * dt / dx

    for i in prange(n_links, nogil=True, schedule="static"):
        link = links[i]

        if fabs(q_at_link[link] * four_dt_over_dx) > h_at_link[link]:
            if q_at_link[link] < 0.0:
                q_at_link[link] = - h_at_link[link] * dx_over_dt * 0.2
            else:
                q_at_link[link] = h_at_link[link] * dx_over_dt * 0.2


@cython.boundscheck(False)
@cython.wraparound(False)
def update_water_depths(
    const cython.floating [:] q_at_node,
    cython.floating [:] h_at_node,
    const cython.floating [:] rainfall_rate_at_node,
    const id_t [:] nodes,
    const double dt,
):
    """
    calculate the change in water depths on all core nodes by finding the
    difference between inputs (rainfall) and the inputs/outputs (flux
    divergence of discharge)
    """
    cdef long n_nodes = len(nodes)
    cdef long i
    cdef long node

    for i in prange(n_nodes, nogil=True, schedule="static"):
        node = nodes[i]

        h_at_node[node] = h_at_node[node] + (
            rainfall_rate_at_node[node] - q_at_node[node]
        ) * dt

        if h_at_node[node] < 0.0:
            h_at_node[node] = 0.0


@cython.boundscheck(False)
@cython.wraparound(False)
def map_sum_of_influx_to_node(
    const cython.floating [:] value_at_link,
    const id_t [:, :] links_at_node,
    const int8_t [:, :] link_dirs_at_node,
    cython.floating [:] out,
):
    """
    calculate the change in water depths on all core nodes by finding the
    difference between inputs (rainfall) and the inputs/outputs (flux
    divergence of discharge)
    """
    cdef long n_nodes = len(links_at_node)
    cdef long links_per_node = links_at_node.shape[1]
    cdef long node
    cdef long col
    cdef double total
    cdef double value

    for node in prange(n_nodes, nogil=True, schedule="static"):
        total = 0.0
        for col in range(links_per_node):
            value = (
                value_at_link[links_at_node[node, col]] * link_dirs_at_node[node, col]
            )
            if value > 0.0:
                total = total + value
        out[node] = total

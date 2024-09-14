cimport cython

from cython.parallel import prange

# https://cython.readthedocs.io/en/stable/src/userguide/fusedtypes.html
ctypedef fused id_t:
    cython.integral
    long long


ctypedef fused integral_out_t:
    cython.integral
    long long


@cython.boundscheck(False)
@cython.wraparound(False)
cpdef void fill_offsets_to_sorted_blocks(
    const id_t [:] array,
    integral_out_t [:] offsets,
) noexcept nogil:
    cdef long n_values = len(array)
    cdef long n_offsets = len(offsets)
    cdef long first_non_negative
    cdef long i
    cdef long j

    for i in range(0, n_values):
        if array[i] >= 0:
            first_non_negative = i
            break
    else:
        first_non_negative = n_values

    offsets[0] = first_non_negative
    for i in range(1, n_offsets):
        offsets[i] = offsets[i - 1]
        for j in range(offsets[i], n_values):
            if array[j] >= i:
                offsets[i] = j
                break
        else:
            offsets[i] = n_values


@cython.boundscheck(False)
@cython.wraparound(False)
cpdef void find_pairs(
    const id_t [:, :] array,
    const integral_out_t [:] offsets,
    const id_t [:, :] pairs,
    id_t [:] where,
) noexcept nogil:
    cdef long n_pairs = len(pairs)
    cdef long pair
    cdef long first
    cdef long second
    cdef long ind

    for pair in prange(n_pairs, nogil=True, schedule="static"):
        first = pairs[pair, 0]
        second = pairs[pair, 1]

        ind = find_pair(array, offsets, first, second)
        if ind == -1:
            ind = find_pair(array, offsets, second, first)

        where[pair] = ind


@cython.boundscheck(False)
@cython.wraparound(False)
cpdef void find_rolling_pairs(
    const id_t [:, :] array,
    const integral_out_t [:] offsets,
    const id_t [:] pairs,
    integral_out_t [:] where,
    const int wraparound,
) noexcept nogil:
    cdef long n_pairs = len(pairs)
    cdef long first
    cdef long second
    cdef long i
    cdef long ind

    for i in range(n_pairs - 1):
        first = pairs[i]
        second = pairs[i + 1]

        ind = find_pair(array, offsets, first, second)
        if ind == -1:
            ind = find_pair(array, offsets, second, first)

        where[i] = ind

    if wraparound:
        first = pairs[n_pairs - 1]
        second = pairs[0]

        ind = find_pair(array, offsets, first, second)
        if ind == -1:
            ind = find_pair(array, offsets, second, first)

        where[n_pairs - 1] = ind


@cython.boundscheck(False)
@cython.wraparound(False)
cpdef void find_rolling_pairs_2d(
    const id_t [:, :] array,
    const integral_out_t [:] offsets,
    const id_t [:, :] pairs,
    integral_out_t [:, :] where,
    const int wraparound,
) noexcept nogil:
    cdef long n_rows = len(pairs)
    cdef long n_cols = pairs.shape[1]
    cdef long length_of_row
    cdef long row

    for row in prange(n_rows, nogil=True, schedule="static"):
        length_of_row = find_first(pairs[row, :], -1)
        if length_of_row == -1:
            length_of_row = n_cols

        if length_of_row > 1:
            find_rolling_pairs(
                array,
                offsets,
                pairs[row, :length_of_row],
                where[row, :length_of_row],
                wraparound,
            )


@cython.boundscheck(False)
@cython.wraparound(False)
cpdef long find_pair(
    const id_t [:, :] pairs,
    const integral_out_t [:] offsets,
    const long a,
    const long b,
) noexcept nogil:
    cdef long n_offsets = len(offsets)
    cdef long start
    cdef long stop
    cdef long ind

    if a + 1 >= n_offsets:
        return -1

    start = offsets[a]
    stop = offsets[a + 1]

    ind = find_first(pairs[start:stop, 1], b)
    if ind >= 0:
        ind += start

    return ind


@cython.boundscheck(False)
@cython.wraparound(False)
cpdef long find_first(
    const id_t [:] array,
    const id_t value,
) noexcept nogil:
    cdef long n_values = len(array)
    cdef long i

    for i in range(n_values):
        if array[i] == value:
            return i
    else:
        return -1

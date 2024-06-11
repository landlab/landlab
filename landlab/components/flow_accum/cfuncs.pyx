import numpy as np

cimport cython
cimport numpy as np

DTYPE_INT = int
ctypedef np.int_t DTYPE_INT_t

DTYPE_FLOAT = np.double
ctypedef np.double_t DTYPE_FLOAT_t


@cython.boundscheck(False)
cpdef _add_to_stack(
    DTYPE_INT_t l,
    DTYPE_INT_t j,
    np.ndarray[DTYPE_INT_t, ndim=1] s,
    np.ndarray[DTYPE_INT_t, ndim=1] delta,
    np.ndarray[DTYPE_INT_t, ndim=1] donors,
):
    """
    Adds node l to the stack and increments the current index (j).
    """
    cdef int m, n, delta_l, delta_lplus1

    s[j] = l
    j += 1
    delta_l = delta[l]
    delta_lplus1 = delta[l+1]

    for n in range(delta_l, delta_lplus1):
        m = donors[n]
        if m != l:
            j = _add_to_stack(m, j, s, delta, donors)

    return j


@cython.boundscheck(False)
cpdef _accumulate_to_n(
    DTYPE_INT_t size,
    DTYPE_INT_t q,
    np.ndarray[DTYPE_INT_t, ndim=1] s,
    np.ndarray[DTYPE_INT_t, ndim=2] r,
    np.ndarray[DTYPE_FLOAT_t, ndim=2] p,
    np.ndarray[DTYPE_FLOAT_t, ndim=1] drainage_area,
    np.ndarray[DTYPE_FLOAT_t, ndim=1] discharge,
):
    """
    Accumulates drainage area and discharge, permitting transmission losses.
    """
    cdef int donor, recvr, i, v
    cdef float accum, proportion

    # Iterate backward through the list, which means we work from upstream to
    # downstream.
    for i in range(size - 1, -1, -1):
        donor = s[i]
        for v in range(q):
            recvr = r[donor, v]
            proportion = p[donor, v]
            if proportion > 0.:
                if donor != recvr:
                    drainage_area[recvr] += proportion*drainage_area[donor]
                    accum = discharge[recvr] + proportion*discharge[donor]
                    if accum < 0.:
                        accum = 0.
                    discharge[recvr] = accum


@cython.boundscheck(False)
cpdef _accumulate_bw(
    DTYPE_INT_t size,
    np.ndarray[DTYPE_INT_t, ndim=1] s,
    np.ndarray[DTYPE_INT_t, ndim=1] r,
    np.ndarray[DTYPE_FLOAT_t, ndim=1] drainage_area,
    np.ndarray[DTYPE_FLOAT_t, ndim=1] discharge,
):
    """
    Accumulates drainage area and discharge, permitting transmission losses.
    """
    cdef int donor, recvr, i
    cdef float accum

    # Iterate backward through the list, which means we work from upstream to
    # downstream.
    for i in range(size - 1, -1, -1):
        donor = s[i]
        recvr = r[donor]
        if donor != recvr:
            drainage_area[recvr] += drainage_area[donor]
            accum = discharge[recvr] + discharge[donor]
            if accum < 0.:
                accum = 0.
            discharge[recvr] = accum


@cython.boundscheck(False)
cpdef _make_donors(
    DTYPE_INT_t size,
    np.ndarray[DTYPE_INT_t, ndim=1] w,
    np.ndarray[DTYPE_INT_t, ndim=1] D,
    np.ndarray[DTYPE_INT_t, ndim=1] delta,
    np.ndarray[DTYPE_INT_t, ndim=1] r,
):
    """Determines number of donors"""
    cdef int ri, i
    for i in range(size):
        ri = r[i]
        D[delta[ri] + w[ri]] = i
        w[ri] += 1


@cython.boundscheck(False)
cpdef _make_donors_to_n(
    DTYPE_INT_t size,
    DTYPE_INT_t q,
    np.ndarray[DTYPE_INT_t, ndim=1] w,
    np.ndarray[DTYPE_INT_t, ndim=1] D,
    np.ndarray[DTYPE_INT_t, ndim=1] delta,
    np.ndarray[DTYPE_INT_t, ndim=2] r,
    np.ndarray[DTYPE_FLOAT_t, ndim=2] p,
):
    """Determines number of donors for route to n"""
    cdef int ri, i, v, ind
    for v in range(q):
        for i in range(size):
            ri = r[i, v]
            if p[i, v] > 0:
                ind = delta[ri] + w[ri]
                D[ind] = i
                w[ri] += 1

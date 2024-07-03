cimport cython


@cython.boundscheck(False)
cpdef _add_to_stack(
    long l,
    long j,
    cython.integral [:] s,
    cython.integral [:] delta,
    cython.integral [:] donors,
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
    long size,
    long q,
    cython.integral [:] s,
    cython.integral [:, :] r,
    cython.floating [:, :] p,
    cython.floating [:] drainage_area,
    cython.floating [:] discharge,
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
    long size,
    cython.integral [:] s,
    cython.integral [:] r,
    cython.floating [:] drainage_area,
    cython.floating [:] discharge,
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
    long size,
    cython.integral [:] w,
    cython.integral [:] D,
    cython.integral [:] delta,
    cython.integral [:] r,
):
    """Determines number of donors"""
    cdef int ri, i
    for i in range(size):
        ri = r[i]
        D[delta[ri] + w[ri]] = i
        w[ri] += 1


@cython.boundscheck(False)
cpdef _make_donors_to_n(
    long size,
    long q,
    cython.integral [:] w,
    cython.integral [:] D,
    cython.integral [:] delta,
    cython.integral [:, :] r,
    cython.floating [:, :] p,
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

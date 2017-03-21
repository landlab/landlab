import numpy as np
cimport numpy as np
cimport cython


DTYPE = np.int
ctypedef np.int_t DTYPE_INT_t


@cython.boundscheck(False)
cpdef _add_to_stack(DTYPE_INT_t l, DTYPE_INT_t j,
                    np.ndarray[DTYPE_INT_t, ndim=1] s,
                    np.ndarray[DTYPE_INT_t, ndim=1] delta,
                    np.ndarray[DTYPE_INT_t, ndim=1] donors):

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

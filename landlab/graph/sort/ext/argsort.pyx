import numpy as np
cimport numpy as np
cimport cython


from libc.stdlib cimport malloc, free
# from libc.stdlib cimport mergesort, qsort


ctypedef np.int_t INT_t
ctypedef np.float_t FLOAT_t

cdef struct Sorter:
    INT_t index
    FLOAT_t value


cdef struct IntSorter:
    INT_t index
    INT_t value


cdef extern from "stdlib.h":
    ctypedef void const_void "const void"
    int qsort(void *, size_t, size_t,
               int(*)(const_void *, const_void *)) nogil
    int	mergesort(void *, size_t, size_t,
                  int (*)(const_void *, const_void *)) nogil


cdef int _compare(const_void *a, const_void *b):
    cdef double v = ((<Sorter*>a)).value - ((<Sorter*>b)).value
    if v < 0:
        return -1
    elif v > 0:
        return 1
    else:
        return 0


cdef int _compare_int(const_void *a, const_void *b):
    cdef int v = ((<IntSorter*>a)).value - ((<IntSorter*>b)).value
    if v < 0:
        return -1
    elif v > 0:
        return 1
    else:
        return 0


cdef void _argsort(double * data, int n_elements, Sorter * order):
    cdef int i

    for i in range(n_elements):
        order[i].index = i
        order[i].value = data[i]

    qsort(<void*> order, n_elements, sizeof(Sorter), _compare)


cdef void _argsort_int(long * data, int n_elements, IntSorter * order):
    cdef int i

    for i in range(n_elements):
        order[i].index = i
        order[i].value = data[i]

    qsort(<void*> order, n_elements, sizeof(IntSorter), _compare_int)
    # mergesort(<void*> order, n_elements, sizeof(IntSorter), _compare_int)


cdef void argsort(double * data, int n_elements, int * out):
    cdef Sorter *sorted_struct = <Sorter*>malloc(n_elements * sizeof(Sorter))

    try:
        _argsort(data, n_elements, sorted_struct)
        for i in range(n_elements):
            out[i] = sorted_struct[i].index
    finally:
        free(sorted_struct)


cdef void argsort_int(long * data, int n_elements, int * out):
    cdef IntSorter *sorted_struct = <IntSorter*>malloc(n_elements *
                                                       sizeof(IntSorter))

    try:
        _argsort_int(data, n_elements, sorted_struct)
        for i in range(n_elements):
            out[i] = sorted_struct[i].index
    finally:
        free(sorted_struct)


cdef int unique_int(long * data, int n_elements, int * out):
    cdef int * index = <int *>malloc(n_elements * sizeof(int))
    cdef int n_unique

    try:
        argsort_int(data, n_elements, index)

        if data[index[0]] != data[index[1]]:
            out[0] = data[index[0]]
        else:
            out[0] = data[index[1]]

        n_unique = 1
        for i in range(1, n_elements):
            if data[index[i]] != data[index[i - 1]]:
                out[n_unique] = data[index[i]]
                n_unique += 1
    finally:
        free(index)

    return n_unique

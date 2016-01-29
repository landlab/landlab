import numpy as np
cimport numpy as np
cimport cython


from libc.stdlib cimport malloc, free, qsort
from libc.math cimport atan2


ctypedef np.int_t INT_t
ctypedef np.float_t FLOAT_t

cdef struct Sorter:
    INT_t index
    FLOAT_t value


cdef extern from "stdlib.h":
    ctypedef void const_void "const void"
    void qsort(void *base, int nmemb, int size,
            int(*compar)(const_void *, const_void *)) nogil


cdef int _compare(const_void *a, const_void *b):
    cdef double v = ((<Sorter*>a)).value - ((<Sorter*>b)).value
    if v < 0:
        return -1
    else:
        return 1


cdef void _argsort(double * data, int n_elements, Sorter * order):
    cdef int i

    for i in range(n_elements):
        order[i].index = i
        order[i].value = data[i]

    qsort(<void*> order, n_elements, sizeof(Sorter), _compare)


cdef void argsort(double * data, int n_elements, int * out):
    cdef Sorter *sorted_struct = <Sorter*>malloc(n_elements * sizeof(Sorter))

    try:
        _argsort(data, n_elements, sorted_struct)
        for i in range(n_elements):
            out[i] = sorted_struct[i].index
    finally:
        free(sorted_struct)

cimport cython

ctypedef fused id_t:
    cython.integral
    long long


cdef void argsort_flt(cython.floating * base, int n_elements, id_t * out) nogil
cdef void argsort_int(long * base, int n_elements, int * out) nogil
cdef int unique_int(long * data, int n_elements, int * out)

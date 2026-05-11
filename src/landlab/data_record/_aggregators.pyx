cimport cython
from cython.parallel cimport prange
from libc.math cimport exp
from libc.math cimport log
from libc.stdint cimport uint8_t
from libc.stdlib cimport free
from libc.stdlib cimport malloc

ctypedef fused id_t:
    cython.integral
    long long


ctypedef fused integral_out_t:
    cython.integral
    long long


ctypedef fused float_or_int:
    cython.integral
    long long
    cython.floating


ctypedef fused float_or_int_weights:
    cython.integral
    long long
    cython.floating


@cython.boundscheck(False)
@cython.wraparound(False)
cpdef void aggregate_items_as_count(
    integral_out_t [:] out,
    const id_t [:] element_of_item,
) noexcept nogil:
    cdef long number_of_elements = len(out)
    cdef long number_of_items = len(element_of_item)
    cdef int item, element

    for element in prange(number_of_elements, nogil=True, schedule="static"):
        out[element] = 0

    for item in range(number_of_items):
        element = element_of_item[item]
        if element >= 0:
            out[element] = out[element] + 1


@cython.boundscheck(False)
@cython.wraparound(False)
cpdef void aggregate_items_as_sum(
    cython.floating [:] out,
    const id_t [:] element_of_item,
    const float_or_int [:] value_of_item,
) noexcept nogil:
    cdef long number_of_elements = len(out)
    cdef long number_of_items = len(element_of_item)
    cdef int item, element

    for element in prange(number_of_elements, nogil=True, schedule="static"):
        out[element] = 0

    for item in range(number_of_items):
        element = element_of_item[item]
        if element >= 0:
            out[element] = out[element] + value_of_item[item]


@cython.boundscheck(False)
@cython.wraparound(False)
cpdef void aggregate_items_as_mean(
    cython.floating [:] out,
    const id_t [:] element_of_item,
    const float_or_int [:] value_of_item,
    const float_or_int_weights [:] weight_of_item,
) noexcept nogil:
    cdef long number_of_elements = len(out)
    cdef long number_of_items = len(element_of_item)
    cdef int item, element
    cdef double * total_weight_at_element = <double *>malloc(
        number_of_elements * sizeof(double)
    )

    try:
        for element in prange(number_of_elements, nogil=True, schedule="static"):
            out[element] = 0.0
            total_weight_at_element[element] = 0.0

        for item in range(number_of_items):
            element = element_of_item[item]
            if element >= 0:
                out[element] = out[element] + value_of_item[item] * weight_of_item[item]
                total_weight_at_element[element] = (
                    total_weight_at_element[element] + weight_of_item[item]
                )

        for element in range(number_of_elements):
            if total_weight_at_element[element] > 0:
                out[element] = out[element] / total_weight_at_element[element]
    finally:
        free(total_weight_at_element)


@cython.boundscheck(False)
@cython.wraparound(False)
@cython.nonecheck(False)
cpdef void aggregate_items_as_gmean(
    cython.floating [::1] out,
    const id_t [::1] element_of_item,
    const float_or_int [::1] value_of_item,
    const float_or_int_weights [::1] weight_of_item,
    const uint8_t[::1] where,
) noexcept nogil:
    cdef Py_ssize_t number_of_elements = len(out)
    cdef Py_ssize_t number_of_items = len(element_of_item)
    cdef Py_ssize_t item
    cdef Py_ssize_t element
    cdef double weight
    cdef double* total_weight = <double*>malloc(
        number_of_elements * sizeof(double)
    )
    cdef double* total_weighted_log = <double*>malloc(
        number_of_elements * sizeof(double)
    )

    try:
        for element in prange(number_of_elements, nogil=True, schedule="static"):
            total_weighted_log[element] = 0.0
            total_weight[element] = 0.0

        for item in range(number_of_items):
            if not where[item]:
                continue

            element = element_of_item[item]
            if element >= 0:
                weight = weight_of_item[item]
                total_weighted_log[element] = (
                    total_weighted_log[element] + weight * log(value_of_item[item])
                )
                total_weight[element] = total_weight[element] + weight

        for element in range(number_of_elements):
            if total_weight[element] > 0.0:
                out[element] = exp(total_weighted_log[element] / total_weight[element])
    finally:
        free(total_weight)
        free(total_weighted_log)

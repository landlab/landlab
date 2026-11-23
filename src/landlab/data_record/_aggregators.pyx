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


ctypedef enum AggregateItemsError:
    OK = 0
    MEMORY_ERROR = 1
    ZERO_TOTAL_WEIGHT_ERROR = 2

AGGREGATE_ITEMS_OK = <uint8_t>OK
AGGREGATE_ITEMS_MEMORY_ERROR = <uint8_t>MEMORY_ERROR
AGGREGATE_ITEMS_ZERO_TOTAL_WEIGHT = <uint8_t>ZERO_TOTAL_WEIGHT_ERROR


@cython.boundscheck(False)
@cython.wraparound(False)
cpdef uint8_t aggregate_items_as_count(
    integral_out_t [::1] out,
    const id_t [::1] element_of_item,
    const uint8_t[::1] where,
) noexcept nogil:
    cdef Py_ssize_t number_of_elements = len(out)
    cdef Py_ssize_t number_of_items = len(element_of_item)
    cdef Py_ssize_t item
    cdef Py_ssize_t element

    for element in prange(number_of_elements, nogil=True, schedule="static"):
        out[element] = 0

    for item in range(number_of_items):
        if not where[item]:
            continue
        element = element_of_item[item]
        out[element] = out[element] + 1

    return OK


@cython.boundscheck(False)
@cython.wraparound(False)
cpdef uint8_t aggregate_items_as_sum(
    cython.floating [::1] out,
    const id_t [::1] element_of_item,
    const float_or_int [::1] value_of_item,
    const uint8_t[::1] where,
) noexcept nogil:
    cdef long number_of_elements = len(out)
    cdef long number_of_items = len(element_of_item)
    cdef int item, element

    for element in prange(number_of_elements, nogil=True, schedule="static"):
        out[element] = 0

    for item in range(number_of_items):
        if not where[item]:
            continue
        element = element_of_item[item]
        out[element] = out[element] + value_of_item[item]

    return OK


@cython.boundscheck(False)
@cython.wraparound(False)
@cython.nonecheck(False)
cpdef uint8_t aggregate_items_as_mean(
    cython.floating [::1] out,
    const id_t [::1] element_of_item,
    const float_or_int [::1] value_of_item,
    const float_or_int_weights [::1] weight_of_item,
    const uint8_t[::1] where,
) noexcept nogil:
    cdef uint8_t status = OK
    cdef Py_ssize_t number_of_elements = len(out)
    cdef Py_ssize_t number_of_items = len(element_of_item)
    cdef Py_ssize_t item
    cdef Py_ssize_t element
    cdef double weight
    cdef double * total_weight = <double *>malloc(number_of_elements * sizeof(double))
    cdef double * total= <double *>malloc(number_of_elements * sizeof(double))
    cdef uint8_t * touched= <uint8_t *>malloc(number_of_elements * sizeof(uint8_t))

    if total_weight == NULL or total == NULL or touched == NULL:
        free(total_weight)
        free(total)
        free(touched)
        return MEMORY_ERROR

    try:
        for element in prange(number_of_elements, nogil=True, schedule="static"):
            total[element] = 0.0
            total_weight[element] = 0.0
            touched[element] = 0

        for item in range(number_of_items):
            if not where[item]:
                continue

            element = element_of_item[item]
            weight = weight_of_item[item]

            total[element] = total[element] + value_of_item[item] * weight
            total_weight[element] = total_weight[element] + weight
            touched[element] = 1

        for element in range(number_of_elements):
            if total_weight[element] > 0:
                out[element] = total[element] / total_weight[element]
            elif touched[element] == 1:
                status = ZERO_TOTAL_WEIGHT_ERROR
    finally:
        free(total)
        free(total_weight)
        free(touched)
    return status


@cython.boundscheck(False)
@cython.wraparound(False)
@cython.nonecheck(False)
cpdef uint8_t aggregate_items_as_gmean(
    cython.floating [::1] out,
    const id_t [::1] element_of_item,
    const float_or_int [::1] value_of_item,
    const float_or_int_weights [::1] weight_of_item,
    const uint8_t[::1] where,
) noexcept nogil:
    cdef uint8_t status = OK
    cdef Py_ssize_t number_of_elements = len(out)
    cdef Py_ssize_t number_of_items = len(element_of_item)
    cdef Py_ssize_t item
    cdef Py_ssize_t element
    cdef double weight
    cdef double* total_weight = <double*>malloc(number_of_elements * sizeof(double))
    cdef double* total = <double*>malloc(number_of_elements * sizeof(double))
    cdef uint8_t * touched = <uint8_t *>malloc(number_of_elements * sizeof(uint8_t))

    if total_weight == NULL or total == NULL or touched == NULL:
        free(total_weight)
        free(total)
        free(touched)
        return MEMORY_ERROR

    try:
        for element in prange(number_of_elements, nogil=True, schedule="static"):
            total[element] = 0.0
            total_weight[element] = 0.0
            touched[element] = 0

        for item in range(number_of_items):
            if not where[item]:
                continue

            element = element_of_item[item]
            weight = weight_of_item[item]

            if weight > 0.0:
                total[element] = total[element] + weight * log(value_of_item[item])
                total_weight[element] = total_weight[element] + weight
            touched[element] = 1

        for element in range(number_of_elements):
            if total_weight[element] > 0.0:
                out[element] = exp(total[element] / total_weight[element])
            elif touched[element]:
                status = ZERO_TOTAL_WEIGHT_ERROR
    finally:
        free(total_weight)
        free(total)
        free(touched)
    return status

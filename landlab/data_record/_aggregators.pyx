import numpy as np

cimport cython
cimport numpy as cnp
from cython.parallel cimport prange
from libc.stdlib cimport free, malloc

ctypedef fused float_or_int:
    cython.integral
    cython.floating


ctypedef fused float_or_int_weights:
    cython.integral
    cython.floating


@cython.boundscheck(False)
@cython.wraparound(False)
def aggregate_items_as_count(
    cython.integral [:] out,
    const long number_of_elements,
    const cython.integral [:] element_of_item,
    const long number_of_items,
):
    cdef int item, element

    for element in prange(number_of_elements, nogil=True, schedule="static"):
        out[element] = 0

    for item in range(number_of_items):
        element = element_of_item[item]
        if element >= 0:
            out[element] = out[element] + 1


@cython.boundscheck(False)
@cython.wraparound(False)
def aggregate_items_as_sum(
    cython.floating [:] out,
    const long number_of_elements,
    const cython.integral [:] element_of_item,
    const long number_of_items,
    const float_or_int [:] value_of_item,
):
    cdef int item, element

    for element in prange(number_of_elements, nogil=True, schedule="static"):
        out[element] = 0

    for item in range(number_of_items):
        element = element_of_item[item]
        if element >= 0:
            out[element] = out[element] + value_of_item[item]


@cython.boundscheck(False)
@cython.wraparound(False)
def aggregate_items_as_mean(
    cython.floating [:] out,
    const long number_of_elements,
    const cython.integral [:] element_of_item,
    const long number_of_items,
    const float_or_int [:] value_of_item,
    const float_or_int_weights [:] weight_of_item,
):
    cdef int item, element
    cdef double * total_weight_at_element = <double *>malloc(number_of_elements * sizeof(double))

    try:
        for element in prange(number_of_elements, nogil=True, schedule="static"):
            out[element] = 0.0
            total_weight_at_element[element] = 0.0

        for item in range(number_of_items):
            element = element_of_item[item]
            if element >= 0:
                out[element] = out[element] + value_of_item[item] * weight_of_item[item]
                total_weight_at_element[element] = total_weight_at_element[element] + weight_of_item[item]

        for element in range(number_of_elements):
            if total_weight_at_element[element] > 0:
                out[element] = out[element] / total_weight_at_element[element]
    finally:
        free(total_weight_at_element)

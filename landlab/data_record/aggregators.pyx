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
def aggregate_items_at_link_count(
    cython.integral [:] out,
    const long number_of_links,
    const cython.integral [:] link_of_item,
    const long number_of_items,
):
    cdef int item, link

    for link in prange(number_of_links, nogil=True, schedule="static"):
        out[link] = 0

    for item in range(number_of_items):
        link = link_of_item[item]
        if link >= 0:
            out[link] = out[link] + 1


@cython.boundscheck(False)
@cython.wraparound(False)
def aggregate_items_at_link_sum(
    cython.floating [:] out,
    const long number_of_links,
    const cython.integral [:] link_of_item,
    const long number_of_items,
    const float_or_int [:] value_of_item,
):
    cdef int item, link

    for link in prange(number_of_links, nogil=True, schedule="static"):
        out[link] = 0

    for item in range(number_of_items):
        link = link_of_item[item]
        if link >= 0:
            out[link] = out[link] + value_of_item[item]


@cython.boundscheck(False)
@cython.wraparound(False)
def aggregate_items_at_link_mean(
    cython.floating [:] out,
    const long number_of_links,
    const cython.integral [:] link_of_item,
    const long number_of_items,
    const float_or_int [:] value_of_item,
    const float_or_int_weights [:] weight_of_item,
):
    cdef int item, link
    cdef double * total_weight_at_link = <double *>malloc(number_of_links * sizeof(double))

    try:
        for link in prange(number_of_links, nogil=True, schedule="static"):
            out[link] = 0.0
            total_weight_at_link[link] = 0.0

        for item in range(number_of_items):
            link = link_of_item[item]
            if link >= 0:
                out[link] = out[link] + value_of_item[item] * weight_of_item[item]
                total_weight_at_link[link] = total_weight_at_link[link] + weight_of_item[item]

        for link in range(number_of_links):
            if total_weight_at_link[link] > 0:
                out[link] = out[link] / total_weight_at_link[link]
    finally:
        free(total_weight_at_link)

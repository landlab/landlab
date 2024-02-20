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
def aggregate_parcels_at_link_count(
    cython.integral [:] out,
    const long number_of_links,
    const cython.integral [:] link_of_parcel,
    const long number_of_parcels,
):
    cdef int parcel, link

    for link in prange(number_of_links, nogil=True, schedule="static"):
        out[link] = 0

    for parcel in range(number_of_parcels):
        link = link_of_parcel[parcel]
        if link >= 0:
            out[link] = out[link] + 1


@cython.boundscheck(False)
@cython.wraparound(False)
def aggregate_parcels_at_link_sum(
    cython.floating [:] out,
    const long number_of_links,
    const cython.integral [:] link_of_parcel,
    const long number_of_parcels,
    const float_or_int [:] value_of_parcel,
):
    cdef int parcel, link

    for link in prange(number_of_links, nogil=True, schedule="static"):
        out[link] = 0

    for parcel in range(number_of_parcels):
        link = link_of_parcel[parcel]
        if link >= 0:
            out[link] = out[link] + value_of_parcel[parcel]


@cython.boundscheck(False)
@cython.wraparound(False)
def aggregate_parcels_at_link_mean(
    cython.floating [:] out,
    const long number_of_links,
    const cython.integral [:] link_of_parcel,
    const long number_of_parcels,
    const float_or_int [:] value_of_parcel,
    const float_or_int_weights [:] weight_of_parcel,
):
    cdef int parcel, link
    cdef double * total_weight_at_link = <double *>malloc(number_of_links * sizeof(double))

    try:
        for link in prange(number_of_links, nogil=True, schedule="static"):
            out[link] = 0.0
            total_weight_at_link[link] = 0.0

        for parcel in range(number_of_parcels):
            link = link_of_parcel[parcel]
            if link >= 0:
                out[link] = out[link] + value_of_parcel[parcel] * weight_of_parcel[parcel]
                total_weight_at_link[link] = total_weight_at_link[link] + weight_of_parcel[parcel]

        for link in range(number_of_links):
            if total_weight_at_link[link] > 0:
                out[link] = out[link] / total_weight_at_link[link]
    finally:
        free(total_weight_at_link)

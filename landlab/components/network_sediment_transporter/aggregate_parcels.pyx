import numpy as np

cimport cython
cimport numpy as cnp
from libc.stdlib cimport free, malloc


@cython.boundscheck(False)
@cython.wraparound(False)
def aggregate_parcels_at_link_sum(
    double [:] out,
    const long number_of_links,
    const double [:] value_of_parcel,
    const long [:] link_of_parcel,
    const long number_of_parcels,
):
    cdef int parcel, link

    for link in range(number_of_links):
        out[link] = 0.0

    for parcel in range(number_of_parcels):
        link = link_of_parcel[parcel]
        if link >= 0:
            out[link] = out[link] + value_of_parcel[parcel]


@cython.boundscheck(False)
@cython.wraparound(False)
def aggregate_parcels_at_link_mean(
    double [:] out,
    const long number_of_links,
    const double [:] value_of_parcel,
    const double [:] weight_of_parcel,
    const long [:] link_of_parcel,
    const long number_of_parcels,
):
    cdef int parcel, link
    cdef double * total_weight_at_link = <double *>malloc(number_of_links * sizeof(double))

    try:
        for link in range(number_of_links):
            out[link] = 0.0

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

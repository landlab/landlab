#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import numpy as np
cimport numpy as np
cimport cython

DTYPE_FLOAT = np.double
ctypedef np.double_t DTYPE_FLOAT_t
DTYPE_INT64 = np.int64
ctypedef np.int64_t DTYPE_INT64_t
DTYPE_UINT8 = np.uint8
ctypedef np.uint8_t DTYPE_UINT8_t


cdef double LARGE_ELEV = 9999999999.0


#cdef find_receiver_and_slope():


def redirect_dirs(
        np.ndarray[DTYPE_INT64_t, ndim=1] liminal_nodes,
        np.ndarray[DTYPE_INT64_t, ndim=2] neighbor_array,
        np.ndarray[DTYPE_INT64_t, ndim=2] link_array,
        np.ndarray[DTYPE_INT64_t, ndim=2] diag_neighbor_array,
        np.ndarray[DTYPE_INT64_t, ndim=2] diag_array,
        np.ndarray[DTYPE_UINT8_t, ndim=1] status_at_node,
        int BC_NODE_IS_CLOSED,
        np.ndarray[DTYPE_FLOAT_t, ndim=1] fill_surface,
        np.ndarray[DTYPE_FLOAT_t, ndim=1] neighbor_lengths,
        np.ndarray[DTYPE_INT64_t, ndim=1] receivers,
        np.ndarray[DTYPE_INT64_t, ndim=1] receiverlinks,
        np.ndarray[DTYPE_FLOAT_t, ndim=1] steepestslopes,
        int has_diags,
        ):

    cdef int liminal
    cdef double min_elev

    for liminal in liminal_nodes:
        min_elev = LARGE_ELEV
        min_link = -1
        # for neighbor_set, link_set in zip(
        #         neighbor_arrays, link_arrays
        # ):
        if True: # ONCE FOR THE BASIC NEIGHBORS
            neighbors = neighbor_array[liminal]
            neighbors_valid = np.not_equal(neighbors, -1)
            closednodes = np.equal(
                status_at_node[neighbors],
                BC_NODE_IS_CLOSED,
            )  # closed BCs can't count
            neighbors_valid[closednodes] = False
            neighbors_to_check = neighbors[neighbors_valid]
            if len(neighbors_to_check) == 0:
                continue
            else:
                min_neighbor_now = np.amin(
                    fill_surface[neighbors_to_check]
                )
                if min_neighbor_now < min_elev:
                    min_elev = min_neighbor_now
                    links_available = link_array[liminal][neighbors_valid]
                    min_link_of_valid = np.argmin(
                        fill_surface[neighbors_to_check]
                    )
                    min_receiver = neighbors_to_check[min_link_of_valid]
                    min_link = links_available[min_link_of_valid]
                    max_grad = (
                        fill_surface[liminal] - min_elev
                    ) / neighbor_lengths[min_link]
                else:
                    pass
        if has_diags: #ONCE FOR THE DIAGONS, IF APPLICABLE
            neighbors = diag_neighbor_array[liminal]
            print('diag nbrs:')
            print(neighbors)
            #SOMEHOW AN ID OF 48 IS SHOWING UP IN A 48-NODE GRID
            neighbors_valid = np.not_equal(neighbors, -1)
            closednodes = np.equal(
                status_at_node[neighbors],
                BC_NODE_IS_CLOSED,
            )  # closed BCs can't count
            neighbors_valid[closednodes] = False
            neighbors_to_check = neighbors[neighbors_valid]
            if len(neighbors_to_check) == 0:
                continue
            else:
                min_neighbor_now = np.amin(
                    fill_surface[neighbors_to_check]
                )
                if min_neighbor_now < min_elev:
                    min_elev = min_neighbor_now
                    links_available = diag_array[liminal][neighbors_valid]
                    min_link_of_valid = np.argmin(
                        fill_surface[neighbors_to_check]
                    )
                    min_receiver = neighbors_to_check[min_link_of_valid]
                    min_link = links_available[min_link_of_valid]
                    max_grad = (
                        fill_surface[liminal] - min_elev
                    ) / neighbor_lengths[min_link]
                else:
                    pass
        assert min_link != -1, neighbors_valid
        # ^link successfully found
        receivers[liminal] = min_receiver
        receiverlinks[liminal] = min_link
        steepestslopes[liminal] = max_grad

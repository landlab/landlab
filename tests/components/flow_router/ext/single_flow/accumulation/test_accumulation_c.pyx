#distutils: language = c++
#distutils: define_macros=NPY_NO_DEPRECATED_API=NPY_1_7_API_VERSION

import numpy as np

cimport cython
cimport numpy as cnp
from libcpp cimport bool

from numpy.testing import assert_array_almost_equal, assert_array_equal

import landlab.components.flow_router.ext.single_flow.accumulation.accumulation as accumulation

from landlab.components.flow_router.ext.single_flow.accumulation cimport (
    accumulation as accumulation_c,
)


def test_add_to_upstream_ordered_nodes():
    # on a grid of 16 nodes
    cdef:
        cnp.int64_t receiver_id = 7, node_idx_in_stack = 0
        cnp.int64_t [:] upstream_ordered_nodes = np.zeros(16, dtype=np.int64)
        cnp.int64_t [:] donors_start_indexes = np.int64([0, 1, 2, 3, 4,
            4, 5, 6, 8, 9])
        cnp.int64_t [:] donors = np.int64([0, 1, 2, 3, 5, 6, 4, 7, 8])

    accumulation_c._add_to_upstream_ordered_nodes(receiver_id, node_idx_in_stack,
                            upstream_ordered_nodes,
                            donors_start_indexes, donors)
    assert receiver_id == 7
    assert node_idx_in_stack == 0
    assert_array_equal(upstream_ordered_nodes, np.int64([7, 4, 0, 0, 0,
            0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]))
    assert_array_equal(donors_start_indexes, np.int64([0, 1, 2, 3, 4,
            4, 5, 6, 8, 9]))
    assert_array_equal(donors, np.int64([0, 1, 2, 3, 5,
            6, 4, 7, 8]))


def test_calc_upstream_order_for_nodes_c():
     # on 25 nodes
    cdef:
        cnp.int64_t[:] base_level_and_closed_nodes = np.int64([0, 1, 2, 4,
            5, 9, 10, 14, 15, 19, 20, 21, 22, 23, 24, 3])
        cnp.int64_t[:] upstream_ordered_nodes = -1 * np.ones(25, dtype=np.int64)
        cnp.int64_t[:] donors_start_indexes = np.int64([0, 1, 3, 5, 6,
            7, 8, 8, 8, 8, 10, 12, 12, 12, 13, 15, 16,
           17, 17, 17, 18, 19, 21, 22, 24, 25])
        cnp.int64_t[:] donors = np.int64([0, 1, 6, 2, 7, 3, 4, 5, 8,
            9, 10, 11, 12, 13, 14, 15, 17,
           19, 20, 16, 21, 22, 18, 23, 24])

    accumulation_c._calc_upstream_order_for_nodes_c(base_level_and_closed_nodes,
            upstream_ordered_nodes, donors_start_indexes, donors)

    assert_array_equal(base_level_and_closed_nodes, np.int64([0, 1, 2, 4,
            5, 9, 10, 14, 15, 19, 20, 21, 22, 23, 24, 3]))
    assert_array_equal(upstream_ordered_nodes, np.int64([0, 1, 6, 2, 7,
        4, 5, 9, 8, 10, 11, 14, 13, 12, 15, 19, 20,
        21, 16, 17, 22, 23, 18, 24, 3]))
    assert_array_equal(donors_start_indexes, np.int64([0, 1, 3, 5, 6,
            7, 8, 8, 8, 8, 10, 12, 12, 12, 13, 15, 16,
           17, 17, 17, 18, 19, 21, 22, 24, 25]))
    assert_array_equal(donors, np.int64([0, 1, 6, 2, 7, 3, 4, 5, 8,
            9, 10, 11, 12, 13, 14, 15, 17,
           19, 20, 16, 21, 22, 18, 23, 24]))

def test_calc_upstream_order_for_nodes():
    # on 25 nodes
    cdef:
        cnp.int64_t[:] base_level_and_closed_nodes = np.int64([0, 1, 2, 4,
            5, 9, 10, 14, 15, 19, 20, 21, 22, 23, 24, 3])
        cnp.int64_t[:] upstream_ordered_nodes = -1 * np.ones(25, dtype=np.int64)
        cnp.int64_t[:] donors_start_indexes = np.int64([0, 1, 3, 5, 6,
            7, 8, 8, 8, 8, 10, 12, 12, 12, 13, 15, 16,
           17, 17, 17, 18, 19, 21, 22, 24, 25])
        cnp.int64_t[:] donors = np.int64([0, 1, 6, 2, 7, 3, 4, 5, 8,
            9, 10, 11, 12, 13, 14, 15, 17,
           19, 20, 16, 21, 22, 18, 23, 24])

    accumulation._calc_upstream_order_for_nodes(base_level_and_closed_nodes,
            upstream_ordered_nodes, donors_start_indexes, donors)

    assert_array_equal(base_level_and_closed_nodes, np.int64([0, 1, 2, 4,
            5, 9, 10, 14, 15, 19, 20, 21, 22, 23, 24, 3]))
    assert_array_equal(upstream_ordered_nodes, np.int64([0, 1, 6, 2, 7,
        4, 5, 9, 8, 10, 11, 14, 13, 12, 15, 19, 20,
        21, 16, 17, 22, 23, 18, 24, 3]))
    assert_array_equal(donors_start_indexes, np.int64([0, 1, 3, 5, 6,
            7, 8, 8, 8, 8, 10, 12, 12, 12, 13, 15, 16,
           17, 17, 17, 18, 19, 21, 22, 24, 25]))
    assert_array_equal(donors, np.int64([0, 1, 6, 2, 7, 3, 4, 5, 8,
            9, 10, 11, 12, 13, 14, 15, 17,
           19, 20, 16, 21, 22, 18, 23, 24]))

def test_calc_drainage_areas():
    # on 25 nodes
    cdef:
        cnp.int64_t [:] downstream_ordered_nodes = np.int64([3, 24, 18, 23, 22,
            17, 16, 21, 20, 19, 15, 12, 13, 14, 11, 10, 8,
            9, 5, 4, 7, 2, 6, 1, 0])
        cnp.int64_t [:] receivers = np.int64([0, 1, 2, 3, 4, 5, 1, 2, 9,
            9, 10, 10, 13, 14, 14, 15, 21,
           16, 23, 19, 20, 21, 22, 23, 24])
        cnp.float64_t [:] drainage_areas = np.float64([0., 0., 0.,
            0., 0., 0., 900., 900., 900.,
            0., 0., 900., 900., 900., 0., 0., 900., 900.,
            900., 0., 0., 0., 0., 0., 0.])

    accumulation._calc_drainage_areas(downstream_ordered_nodes, receivers,
                                                drainage_areas)
    assert_array_equal(downstream_ordered_nodes, np.int64([3, 24, 18, 23, 22,
            17, 16, 21, 20, 19, 15, 12, 13, 14, 11, 10, 8,
            9, 5, 4, 7, 2, 6, 1, 0]))
    assert_array_equal(receivers, np.int64([0, 1, 2, 3, 4, 5, 1, 2, 9,
            9, 10, 10, 13, 14, 14, 15, 21,
           16, 23, 19, 20, 21, 22, 23, 24]))
    assert_array_almost_equal(drainage_areas, np.float64([0., 900., 900.,
            0., 0., 0., 900., 900.,
            900., 900., 900., 900., 900., 1800., 1800., 0.,
            1800., 900., 900., 0., 0., 1800., 0., 900.,
            0.]))

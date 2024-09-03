# distutils: define_macros=NPY_NO_DEPRECATED_API=NPY_1_7_API_VERSION
# distutils: extra_compile_args = -std=c++11
# distutils: extra_link_args = -std=c++11

# NB: apparently not possible to add language: C++ in this file
# because of the extracompile -std=c++11 (necessary to understand the
# priorityqueue template. To be compiled in C++,
# must add a .pxd file with the instruction # distutils: language = c++

import numpy as np

cimport numpy as cnp
from libcpp cimport bool
from libcpp.pair cimport pair

from numpy.testing import assert_array_almost_equal
from numpy.testing import assert_array_equal

from landlab.components.flow_router.ext.single_flow.priority_routing cimport (
    breach as breach_c,
)

import landlab.components.flow_router.ext.single_flow.priority_routing.breach as breach


cdef extern from "_priority_queue.hpp" nogil:
    cdef cppclass _priority_queue:
        _priority_queue(...) except +
        void push(pair[cnp.int_t, cnp.float_t])
        pair[cnp.int_t, cnp.float_t] top() except +
        void pop()
        bool empty()
        cnp.int_t size()


def test_priority_queue():
    cdef:
        pair[cnp.int_t, cnp.float_t] a = pair[cnp.int_t, cnp.float_t](0, 1045.3)
        pair[cnp.int_t, cnp.float_t] b = pair[cnp.int_t, cnp.float_t](1, 536.3)
        pair[cnp.int_t, cnp.float_t] c = pair[cnp.int_t, cnp.float_t](2, 2034.12)
        _priority_queue to_do = _priority_queue(breach_c._compare_second)

    assert to_do.empty() is True
    to_do.push(a)
    to_do.push(b)
    assert to_do.empty() is False
    assert to_do.top() == b
    to_do.push(c)
    to_do.pop()
    assert to_do.top() == a


def test_init_flow_direction_queues():
    # on a grid of 7 nodes
    cdef:
        cnp.int_t nodes_n = 7
        cnp.int_t [:] base_level_nodes = np.array([0, 2])
        cnp.int_t [:] closed_nodes = np.array([5])
        cnp.float_t [:] z = np.array([4.5, 3.2, 6.7, 13.2, 5.6, 100.3, 45.32])
        _priority_queue to_do = _priority_queue(breach_c._compare_second)
        cnp.int_t [:] receivers = -1 * np.ones(nodes_n, dtype=int)
        cnp.int_t [:] outlet_nodes = -1 * np.ones(nodes_n, dtype=int)
        cnp.int_t [:] done = np.zeros(nodes_n, dtype=int)
        cnp.int_t done_n = 0
    breach_c._init_flow_direction_queues(
        base_level_nodes, closed_nodes, z, to_do, receivers, outlet_nodes, done, &done_n
    )
    assert nodes_n == 7
    assert_array_equal(base_level_nodes, np.array([0, 2]))
    assert_array_equal(closed_nodes, np.array([5]))
    assert_array_almost_equal(z, np.array([4.5, 3.2, 6.7, 13.2, 5.6, 100.3, 45.32]))
    assert_array_equal(receivers, np.array([0, -1, 2, -1, -1, 5, -1]))
    assert_array_equal(done, np.array([1, 0, 1, 0, 0, 1, 0]))
    assert done_n == 3


def test_set_flooded_and_outlet():
    # on a grid of 5 nodes
    cdef:
        cnp.int_t donor_id = 2, bad_index = -1, flooded_status = 3
        cnp.float_t min_elevation_relative_diff = 1e-2
        cnp.float_t [:] z = np.array([0.1, 67.1, 42.1, 70.3, 34.5])
        cnp.int_t [:] receivers = np.array([0, 0, 1, -1, -1])
        cnp.int_t [:] outlet_nodes = np.array([0, 0, -1, -1, -1])
        cnp.int_t [:] depression_outlet_nodes = np.array([-1, -1, -1, -1, -1])
        cnp.int_t [:] flooded_nodes = np.zeros(5, dtype=int)
        cnp.float_t [:] depression_depths = np.zeros(5, dtype=float)
        cnp.float_t [:] depression_free_elevations = (
            np.array([0.1, 67.1, 42.1, 70.3, 34.5])
        )

    breach_c._set_flooded_and_outlet(
        donor_id,
        z,
        receivers,
        outlet_nodes,
        depression_outlet_nodes,
        flooded_nodes,
        depression_depths,
        depression_free_elevations,
        flooded_status,
        bad_index,
        min_elevation_relative_diff,
    )

    assert donor_id == 2
    assert_array_almost_equal(z, np.array([0.1, 67.1, 42.1, 70.3, 34.5]))
    assert_array_equal(receivers, np.array([0, 0, 1, -1, -1]))
    assert_array_equal(outlet_nodes, np.array([0, 0, 0, -1, -1]))
    assert_array_equal(depression_outlet_nodes, np.array([-1, -1, 1, -1, -1]))
    assert_array_equal(flooded_nodes, np.array([0, 0, 3, 0, 0]))
    assert_array_almost_equal(depression_depths, np.array([0, 0, 67.1 - 42.1, 0, 0]))
    assert_array_almost_equal(
        depression_free_elevations, np.array([0.1, 67.1, 67.771, 70.3, 34.5])
    )
    assert flooded_status == 3
    assert bad_index == -1

    donor_id = 3
    receivers = np.array([0, 0, 1, 2, -1])

    breach_c._set_flooded_and_outlet(
        donor_id,
        z,
        receivers,
        outlet_nodes,
        depression_outlet_nodes,
        flooded_nodes,
        depression_depths,
        depression_free_elevations,
        flooded_status,
        bad_index,
        min_elevation_relative_diff,
    )
    assert donor_id == 3
    assert_array_almost_equal(z, np.array([0.1, 67.1, 42.1, 70.3, 34.5]))
    assert_array_equal(receivers, np.array([0, 0, 1, 2, -1]))
    assert_array_equal(outlet_nodes, np.array([0, 0, 0, 0, -1]))
    assert_array_equal(depression_outlet_nodes, np.array([-1, -1, 1, -1, -1]))
    assert_array_equal(flooded_nodes, np.array([0, 0, 3, 0, 0]))
    assert_array_almost_equal(depression_depths, np.array([0, 0, 67.1 - 42.1, 0, 0]))
    assert_array_almost_equal(
        depression_free_elevations, np.array([0.1, 67.1, 67.771, 70.3, 34.5])
    )
    assert flooded_status == 3
    assert bad_index == -1


def test_set_receiver():
    # on a grid of 6 nodes
    cdef:
        cnp.int_t donor_id = 3, receiver_id = 4, done_n = 1
        cnp.int_t [:] receivers = np.array([0, -1, -1, -1, -1, -1])
        cnp.int_t [:] done = np.array([1, 0, 0, 0, 0, 0])
    breach_c._set_receiver(donor_id, receiver_id, receivers, done, &done_n)
    assert donor_id == 3
    assert receiver_id == 4
    assert done_n == 2
    assert_array_equal(receivers, np.array([0, -1, -1, 4, -1, -1]))
    assert_array_equal(done, np.array([1, 0, 0, 1, 0, 0]))


def test_set_donor_properties():
    # on a grid of 9 nodes
    cdef:
        cnp.int_t donor_id = 5, receiver_id = 7
        cnp.int_t [:] _sorted_pseudo_heads = np.array(
            [
                0, 0, 1, 1, 1, 2, 2, 3, 3,
                3, 4, 4, 4, 4, 5, 5, 5, 6, 6, 7, 7, 7, 8, 8
            ]
        )
        cnp.int_t [:] sorted_pseudo_tails = np.array(
            [
                3, 1, 0, 4, 2, 1, 5, 0, 6,
                4, 7, 3, 1, 5, 4, 2, 8, 7, 3, 8, 4, 6, 5, 7
            ]
        )
        cnp.int_t [:, :] head_start_end_indexes = np.array(
            [
                [0, 2, 5, 7, 10, 14, 17, 19, 22],
                [1, 4, 6, 9, 13, 16, 18, 21, 23]
            ]
        )
        cnp.int_t [:] sorted_dupli_links = np.array(
            [
                2, 0, 0, 3, 1, 1, 4,
                2, 7, 5, 8, 5, 3, 6, 6, 4, 9,
                10, 7, 11, 8, 10, 9, 11
            ]
        )
        cnp.float_t [:] sorted_dupli_gradients = np.array(
            [
                0.03879335, 0.04387396,
                0.04387396, 0.08236696, 0.12232775,
                0.12232775, 0.02936549, 0.03879335, 0.28953386, 0.16503428,
                0.18134194, 0.16503428, 0.08236696, 0.06932627, 0.06932627,
                0.02936549, 0.00526324, 0.30584151, 0.28953386, 0.24540497,
                0.18134194, 0.30584151, 0.00526324, 0.24540497
            ]
        )
        cnp.float_t [:] z = np.array(
            [
                2.29047865, 3.60669759, 7.27652998, 1.12667805,
                6.0777065, 8.15749462, 9.81269383, 0.63744841, 7.99959748
            ]
        )
        cnp.float_t [:] steepest_slopes = np.zeros(9, dtype=float)
        cnp.int_t [:] links_to_receivers = -1 * np.ones(9, dtype=int)

    breach_c._set_donor_properties(
        donor_id,
        receiver_id,
        sorted_pseudo_tails,
        head_start_end_indexes,
        sorted_dupli_links,
        sorted_dupli_gradients,
        z,
        steepest_slopes,
        links_to_receivers,
    )
    assert donor_id == 5
    assert receiver_id == 7
    assert_array_equal(
        sorted_pseudo_tails,
        np.array(
            [
                3, 1, 0, 4, 2, 1, 5, 0, 6,
                4, 7, 3, 1, 5, 4, 2, 8, 7, 3, 8, 4, 6, 5, 7
            ]
        )
    )
    assert_array_equal(
        head_start_end_indexes,
        np.array(
            [
                [0, 2, 5, 7, 10, 14, 17, 19, 22],
                [1, 4, 6, 9, 13, 16, 18, 21, 23]
            ]
        )
    )
    assert_array_equal(
        sorted_dupli_links,
        np.array(
            [
                2, 0, 0, 3, 1, 1, 4,
                2, 7, 5, 8, 5, 3, 6, 6, 4, 9,
                10, 7, 11, 8, 10, 9, 11
            ]
        )
    )
    assert_array_almost_equal(
        sorted_dupli_gradients,
        np.array(
            [
                0.03879335,
                0.04387396,
                0.04387396, 0.08236696, 0.12232775,
                0.12232775, 0.02936549, 0.03879335, 0.28953386, 0.16503428,
                0.18134194, 0.16503428, 0.08236696, 0.06932627, 0.06932627,
                0.02936549, 0.00526324, 0.30584151, 0.28953386, 0.24540497,
                0.18134194, 0.30584151, 0.00526324, 0.24540497
            ]
        )
    )
    assert_array_almost_equal(
        z,
        np.array(
            [
                2.29047865, 3.60669759, 7.27652998, 1.12667805,
                6.0777065, 8.15749462, 9.81269383, 0.63744841, 7.99959748
            ]
        )
    )
    assert_array_almost_equal(
        steepest_slopes,
        np.array(
            [
                0., 0.,
                0., 0., 0., 0.245405,
                0., 0., 0.
            ]
        )
    )
    assert_array_equal(
        links_to_receivers, np.array([-1, -1, -1, -1, -1, 11, -1, -1, -1])
    )

#######################################################################################


def test_direct_flow_c():
    # Grid of 25 nodes
    cdef:
        cnp.int_t nodes_n = 25, flooded_status = 3, bad_index = -1
        cnp.float_t min_elevation_relative_diff = 1e-2
        cnp.int_t neighbors_max_number = 50
        cnp.int_t[:] base_level_nodes = np.array(
            [0, 1, 2, 4, 5, 9, 10, 14, 15, 19, 20, 21, 22, 23, 24]
        )
        cnp.int_t[:] base_level_nodes_0 = np.copy(base_level_nodes)
        cnp.int_t[:] closed_nodes = np.array([3])
        cnp.int_t[:] closed_nodes_0 = np.copy(closed_nodes)
        cnp.int_t[:] _sorted_pseudo_heads = np.array(
            [
                0, 0, 1, 1, 1,
                2, 2, 2, 3, 3, 3, 4, 4, 5, 5, 5, 6,
                6, 6, 6, 7, 7, 7, 7, 8, 8, 8, 8, 9, 9, 9, 10, 10, 10,
                11, 11, 11, 11, 12, 12, 12, 12, 13, 13, 13, 13, 14, 14, 14, 15, 15,
                15, 16, 16, 16, 16, 17, 17, 17, 17, 18, 18, 18, 18, 19, 19, 19, 20,
                20, 21, 21, 21, 22, 22, 22, 23, 23, 23, 24, 24
            ]
        )  # unused, for info
        cnp.int_t[:] sorted_pseudo_tails = np.array(
            [
                5, 1, 6, 2, 0, 3,
                7, 1, 8, 2, 4, 9, 3, 6, 10, 0, 7,
                11, 1, 5, 12, 8, 6, 2, 7, 9, 3, 13, 8, 4, 14, 15, 11, 5,
                16, 6, 12, 10, 13, 17, 7, 11, 14, 18, 8, 12, 19, 9, 13, 16, 20,
                10, 21, 15, 17, 11, 12, 22, 16, 18, 19, 23, 17, 13, 24, 18, 14, 21,
                15, 20, 22, 16, 17, 23, 21, 24, 18, 22, 19, 23
            ]
        )
        cnp.int_t[:] sorted_pseudo_tails_0 = np.copy(sorted_pseudo_tails)
        cnp.float_t[:] sorted_dupli_gradients = np.array(
            [
                0.1955672,
                0.04387396, 0.20686654, 0.12232775, 0.04387396,
                0.20499506, 0.22130272, 0.12232775, 0.22909731, 0.20499506,
                0.16503428, 0.17027403, 0.16503428, 0.05517331, 0.1382939,
                0.1955672,  0.30584151, 0.10637057, 0.20686654, 0.05517331,
                0.30061369, 0.24540497, 0.30584151, 0.22130272, 0.24540497,
                0.23433706, 0.22909731, 0.10635323, 0.23433706, 0.17027403,
                0.05572124, 0.06218312, 0.08709664, 0.1382939,  0.20119303,
                0.10637057, 0.10114275, 0.08709664, 0.16156195, 0.12329631,
                0.30061369, 0.10114275, 0.07226259, 0.01279225, 0.10635323,
                0.16156195, 0.08976513, 0.05572124, 0.07226259, 0.17627952,
                0.00161023, 0.06218312, 0.1169648,  0.17627952, 0.17903947,
                0.20119303, 0.12329631, 0.0359346,  0.17903947, 0.05105789,
                0.0302948,  0.13452944, 0.05105789, 0.01279225, 0.14332401,
                0.0302948,  0.08976513, 0.06092495, 0.00161023, 0.06092495,
                0.09800927, 0.1169648,  0.0359346,  0.22152193, 0.09800927,
                0.02150023, 0.13452944, 0.22152193, 0.14332401, 0.02150023
            ]
        )
        cnp.float_t[:] sorted_dupli_gradients_0 = np.copy(sorted_dupli_gradients)
        cnp.int_t[:] sorted_dupli_links = np.array(
            [
                4, 0, 5, 1, 0,
                2, 6, 1, 7, 2, 3, 8, 3, 9, 13, 4, 10,
                14, 5, 9, 15, 11, 10, 6, 11, 12, 7, 16, 12, 8, 17, 22, 18, 13,
                23, 14, 19, 18, 20, 24, 15, 19, 21, 25, 16, 20, 26, 17, 21, 27, 31,
                22, 32, 27, 28, 23, 24, 33, 28, 29, 30, 34, 29, 25, 35, 30, 26, 36,
                31, 36, 37, 32, 33, 38, 37, 39, 34, 38, 35, 39
            ]
        )
        cnp.int_t[:] sorted_dupli_links_0 = np.copy(sorted_dupli_links)
        cnp.int_t[:, :] head_start_end_indexes = np.array(
            [
                [
                    0, 2, 5,
                    8, 11, 13, 16, 20, 24, 28, 31, 34, 38, 42, 46, 49, 52,
                    56, 60, 64, 67, 69, 72, 75, 78
                ],
                [
                    1, 4, 7, 10, 12, 15, 19, 23, 27, 30, 33, 37, 41, 45, 48, 51,
                    55, 59, 63, 66, 68, 71, 74, 77, 79
                ]
            ]
        )
        cnp.int_t[:, :] head_start_end_indexes_0 = np.copy(head_start_end_indexes)
        cnp.float_t[:] depression_depths = np.zeros(25, dtype=float)
        cnp.int_t[:] outlet_nodes = -1 * np.ones(25, dtype=int)
        cnp.int_t[:] depression_outlet_nodes = -1 * np.ones(25, dtype=int)
        cnp.int_t[:] flooded_nodes = np.zeros(25, dtype=int)
        cnp.int_t[:] links_to_receivers = -1 * np.ones(25, dtype=int)
        cnp.int_t[:] receivers = -1 * np.ones(25, dtype=int)
        cnp.float_t[:] steepest_slopes = np.zeros(25, dtype=float)
        cnp.float_t[:] z = np.array(
            [
                2.29047865, 3.60669759, 7.27652998,
                1.12667805, 6.0777065,
                8.15749462, 9.81269383, 0.63744841, 7.99959748, 0.96948555,
                4.00867748, 6.62157669, 9.65585909, 4.80900059, 2.64112287,
                5.8741712,  0.5857857,  5.95696968, 4.42523296, 5.33407687,
                5.92247815, 4.09472964, 7.03500768, 0.38934984, 1.03435662
            ]
        )
        cnp.float_t[:] depression_free_elevations = z.copy()
        cnp.float_t[:] z_0 = np.copy(z)

    breach_c._direct_flow_c(
        nodes_n,
        base_level_nodes,
        closed_nodes,
        sorted_pseudo_tails,
        sorted_dupli_gradients,
        sorted_dupli_links,
        head_start_end_indexes,
        outlet_nodes,
        depression_outlet_nodes,
        flooded_nodes,
        depression_depths,
        depression_free_elevations,
        links_to_receivers,
        receivers,
        steepest_slopes,
        z,
        flooded_status,
        bad_index,
        neighbors_max_number,
        min_elevation_relative_diff,
    )
    assert nodes_n == 25
    assert flooded_status == 3
    assert bad_index == -1
    assert neighbors_max_number == 50
    assert_array_equal(base_level_nodes, base_level_nodes_0)
    assert_array_equal(closed_nodes, closed_nodes_0)
    assert_array_equal(sorted_pseudo_tails, sorted_pseudo_tails_0)
    assert_array_almost_equal(sorted_dupli_gradients, sorted_dupli_gradients_0)
    assert_array_equal(sorted_dupli_links, sorted_dupli_links_0)
    assert_array_equal(head_start_end_indexes, head_start_end_indexes_0)
    assert_array_almost_equal(
        depression_depths,
        np.array(
            [
                0., 0., 0., 0., 0.,
                0., 0., 6.63908157, 0., 0.,
                0., 0., 0., 0., 0.,
                0., 3.50894394, 0., 0., 0.,
                0., 0., 0., 0., 0.
            ]
        )
    )
    assert_array_equal(
        outlet_nodes,
        np.array(
            [
                0, 1, 2, 3, 4, 5, 1, 2, 9, 9, 10, 10, 14,
                14, 14, 15, 21, 21, 23, 19, 20, 21, 22, 23, 24
            ]
        )
    )
    assert_array_equal(
        depression_outlet_nodes,
        np.array(
            [
                -1, -1, -1, -1, -1, -1, -1, 2, -1, -1, -1, -1, -1,
                -1, -1, -1, 21, -1, -1, -1, -1, -1, -1, -1, -1
            ]
        )
    )
    assert_array_equal(
        flooded_nodes,
        np.array(
            [
                0, 0, 0, 0, 0, 0, 0, 3, 0, 0, 0, 0, 0, 0, 0, 0, 3,
                0, 0, 0, 0, 0, 0, 0, 0
            ]
        )
    )
    assert_array_equal(
        links_to_receivers,
        np.array(
            [
                -1, -1, -1, -1, -1, -1, 5, 6, 12, -1, -1, 18, 20,
                21, -1, -1, 32, 28, 34, -1, -1, -1, -1, -1, -1
            ]
        )
    )
    assert_array_equal(
        receivers,
        np.array(
            [
                0, 1, 2, 3, 4, 5, 1, 2, 9, 9, 10, 10,
                13, 14, 14, 15, 21, 16, 23, 19, 20, 21, 22, 23, 24
            ]
        )
    )
    assert_array_almost_equal(
        steepest_slopes,
        np.array(
            [
                0., 0., 0., 0., 0.,
                0., 0.20686654, 0., 0.23433706, 0.,
                0., 0.08709664, 0.16156195, 0.07226259, 0.,
                0., 0., 0.17903947, 0.13452944, 0.,
                0., 0., 0., 0., 0.
            ]
        )
    )
    assert_array_almost_equal(z, z_0)


def test_direct_flow():
    # Grid of 25 nodes
    cdef:
        cnp.int_t nodes_n = 25, flooded_status = 3, bad_index = -1
        cnp.int_t neighbors_max_number = 50
        cnp.float_t min_elevation_relative_diff = 1e-2
        cnp.int_t[:] base_level_nodes = np.array(
            [0, 1, 2, 4, 5, 9, 10, 14, 15, 19, 20, 21, 22, 23, 24]
        )
        cnp.int_t[:] base_level_nodes_0 = np.copy(base_level_nodes)
        cnp.int_t[:] closed_nodes = np.array([3])
        cnp.int_t[:] closed_nodes_0 = np.copy(closed_nodes)
        cnp.int_t[:] _sorted_pseudo_heads = np.array(
            [
                0, 0, 1, 1, 1,
                2, 2, 2, 3, 3, 3, 4, 4, 5, 5, 5, 6,
                6, 6, 6, 7, 7, 7, 7, 8, 8, 8, 8, 9, 9, 9, 10, 10, 10,
                11, 11, 11, 11, 12, 12, 12, 12, 13, 13, 13, 13, 14, 14, 14, 15, 15,
                15, 16, 16, 16, 16, 17, 17, 17, 17, 18, 18, 18, 18, 19, 19, 19, 20,
                20, 21, 21, 21, 22, 22, 22, 23, 23, 23, 24, 24
            ]
        )  # unused, for info
        cnp.int_t[:] sorted_pseudo_tails = np.array(
            [
                5, 1, 6, 2, 0, 3,
                7, 1, 8, 2, 4, 9, 3, 6, 10, 0, 7,
                11, 1, 5, 12, 8, 6, 2, 7, 9, 3, 13, 8, 4, 14, 15, 11, 5,
                16, 6, 12, 10, 13, 17, 7, 11, 14, 18, 8, 12, 19, 9, 13, 16, 20,
                10, 21, 15, 17, 11, 12, 22, 16, 18, 19, 23, 17, 13, 24, 18, 14, 21,
                15, 20, 22, 16, 17, 23, 21, 24, 18, 22, 19, 23
            ]
        )
        cnp.int_t[:] sorted_pseudo_tails_0 = np.copy(sorted_pseudo_tails)
        cnp.float_t[:] sorted_dupli_gradients = np.array(
            [
                0.1955672,
                0.04387396, 0.20686654, 0.12232775, 0.04387396,
                0.20499506, 0.22130272, 0.12232775, 0.22909731, 0.20499506,
                0.16503428, 0.17027403, 0.16503428, 0.05517331, 0.1382939,
                0.1955672,  0.30584151, 0.10637057, 0.20686654, 0.05517331,
                0.30061369, 0.24540497, 0.30584151, 0.22130272, 0.24540497,
                0.23433706, 0.22909731, 0.10635323, 0.23433706, 0.17027403,
                0.05572124, 0.06218312, 0.08709664, 0.1382939,  0.20119303,
                0.10637057, 0.10114275, 0.08709664, 0.16156195, 0.12329631,
                0.30061369, 0.10114275, 0.07226259, 0.01279225, 0.10635323,
                0.16156195, 0.08976513, 0.05572124, 0.07226259, 0.17627952,
                0.00161023, 0.06218312, 0.1169648,  0.17627952, 0.17903947,
                0.20119303, 0.12329631, 0.0359346,  0.17903947, 0.05105789,
                0.0302948,  0.13452944, 0.05105789, 0.01279225, 0.14332401,
                0.0302948,  0.08976513, 0.06092495, 0.00161023, 0.06092495,
                0.09800927, 0.1169648,  0.0359346,  0.22152193, 0.09800927,
                0.02150023, 0.13452944, 0.22152193, 0.14332401, 0.02150023
            ]
        )
        cnp.float_t[:] sorted_dupli_gradients_0 = np.copy(sorted_dupli_gradients)
        cnp.int_t[:] sorted_dupli_links = np.array(
            [
                4, 0, 5, 1, 0,
                2, 6, 1, 7, 2, 3, 8, 3, 9, 13, 4, 10,
                14, 5, 9, 15, 11, 10, 6, 11, 12, 7, 16, 12, 8, 17, 22, 18, 13,
                23, 14, 19, 18, 20, 24, 15, 19, 21, 25, 16, 20, 26, 17, 21, 27, 31,
                22, 32, 27, 28, 23, 24, 33, 28, 29, 30, 34, 29, 25, 35, 30, 26, 36,
                31, 36, 37, 32, 33, 38, 37, 39, 34, 38, 35, 39
            ]
        )
        cnp.int_t[:] sorted_dupli_links_0 = np.copy(sorted_dupli_links)
        cnp.int_t[:, :] head_start_end_indexes = np.array(
            [
                [
                    0, 2, 5,
                    8, 11, 13, 16, 20, 24, 28, 31, 34, 38, 42, 46, 49, 52,
                    56, 60, 64, 67, 69, 72, 75, 78
                ],
                [
                    1, 4, 7, 10, 12, 15, 19, 23, 27, 30, 33, 37, 41, 45, 48, 51,
                    55, 59, 63, 66, 68, 71, 74, 77, 79
                ]
            ]
        )
        cnp.int_t[:, :] head_start_end_indexes_0 = np.copy(head_start_end_indexes)
        cnp.float_t[:] depression_depths = np.zeros(25, dtype=float)
        cnp.int_t[:] outlet_nodes = -1 * np.ones(25, dtype=int)
        cnp.int_t[:] depression_outlet_nodes = -1 * np.ones(25, dtype=int)
        cnp.int_t[:] flooded_nodes = np.zeros(25, dtype=int)
        cnp.int_t[:] links_to_receivers = -1 * np.ones(25, dtype=int)
        cnp.int_t[:] receivers = -1 * np.ones(25, dtype=int)
        cnp.float_t[:] steepest_slopes = np.zeros(25, dtype=float)
        cnp.float_t[:] z = np.array(
            [
                2.29047865, 3.60669759, 7.27652998,
                1.12667805, 6.0777065,
                8.15749462, 9.81269383, 0.63744841, 7.99959748, 0.96948555,
                4.00867748, 6.62157669, 9.65585909, 4.80900059, 2.64112287,
                5.8741712,  0.5857857,  5.95696968, 4.42523296, 5.33407687,
                5.92247815, 4.09472964, 7.03500768, 0.38934984, 1.03435662
            ]
        )
        cnp.float_t[:] depression_free_elevations = z.copy()
        cnp.float_t[:] z_0 = np.copy(z)

    breach._direct_flow(
        nodes_n,
        base_level_nodes,
        closed_nodes,
        sorted_pseudo_tails,
        sorted_dupli_gradients,
        sorted_dupli_links,
        head_start_end_indexes,
        outlet_nodes,
        depression_outlet_nodes,
        flooded_nodes,
        depression_depths,
        depression_free_elevations,
        links_to_receivers,
        receivers,
        steepest_slopes,
        z,
        flooded_status,
        bad_index,
        neighbors_max_number,
        min_elevation_relative_diff,
    )
    assert nodes_n == 25
    assert flooded_status == 3
    assert bad_index == -1
    assert neighbors_max_number == 50
    assert_array_equal(base_level_nodes, base_level_nodes_0)
    assert_array_equal(closed_nodes, closed_nodes_0)
    assert_array_equal(sorted_pseudo_tails, sorted_pseudo_tails_0)
    assert_array_almost_equal(sorted_dupli_gradients, sorted_dupli_gradients_0)
    assert_array_equal(sorted_dupli_links, sorted_dupli_links_0)
    assert_array_equal(head_start_end_indexes, head_start_end_indexes_0)
    assert_array_almost_equal(
        depression_depths,
        np.array(
            [
                0., 0., 0., 0., 0.,
                0., 0., 6.63908157, 0., 0.,
                0., 0., 0., 0., 0.,
                0., 3.50894394, 0., 0., 0.,
                0., 0., 0., 0., 0.
            ]
        )
    )
    assert_array_equal(
        outlet_nodes,
        np.array(
            [
                0, 1, 2, 3, 4, 5, 1, 2, 9, 9, 10, 10, 14,
                14, 14, 15, 21, 21, 23, 19, 20, 21, 22, 23, 24
            ]
        )
    )
    assert_array_equal(
        depression_outlet_nodes,
        np.array(
            [
                -1, -1, -1, -1, -1, -1, -1, 2, -1, -1, -1, -1, -1,
                -1, -1, -1, 21, -1, -1, -1, -1, -1, -1, -1, -1
            ]
        )
    )
    assert_array_equal(
        flooded_nodes,
        np.array(
            [
                0, 0, 0, 0, 0, 0, 0, 3, 0, 0, 0, 0, 0, 0, 0, 0, 3,
                0, 0, 0, 0, 0, 0, 0, 0
            ]
        )
    )
    assert_array_equal(
        links_to_receivers,
        np.array(
            [
                -1, -1, -1, -1, -1, -1, 5, 6, 12, -1, -1, 18, 20,
                21, -1, -1, 32, 28, 34, -1, -1, -1, -1, -1, -1
            ]
        )
    )
    assert_array_equal(
        receivers,
        np.array(
            [
                0, 1, 2, 3, 4, 5, 1, 2, 9, 9, 10, 10,
                13, 14, 14, 15, 21, 16, 23, 19, 20, 21, 22, 23, 24
            ]
        )
    )
    assert_array_almost_equal(
        steepest_slopes,
        np.array(
            [
                0., 0., 0., 0., 0.,
                0., 0.20686654, 0., 0.23433706, 0.,
                0., 0.08709664, 0.16156195, 0.07226259, 0.,
                0., 0., 0.17903947, 0.13452944, 0.,
                0., 0., 0., 0., 0.
            ]
        )
    )
    assert_array_almost_equal(z, z_0)

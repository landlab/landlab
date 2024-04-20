#distutils: language = c++
#distutils: define_macros=NPY_NO_DEPRECATED_API=NPY_1_7_API_VERSION

import numpy as np

cimport cython
cimport numpy as cnp

from numpy.testing import assert_array_equal

import landlab.components.flow_router.ext.single_flow.priority_routing.init_tools as init_tools


def test_get_start_end_indexes_in_sorted_array():
    cdef:
        cnp.int_t [:] sorted_array = np.array([0, 0, 0, 1, 2, 2, 3, 3, 3, 3])
        cnp.int_t n_values = 4, max_value = 9
        cnp.int_t [:, :] idx2 = np.empty((2, n_values), dtype=int)

    idx2 = init_tools._get_start_end_indexes_in_sorted_array(sorted_array,
                                                             n_values, max_value)
    assert_array_equal(idx2, np.array([[0, 3, 4, 6], [2, 3, 5, 9]]))

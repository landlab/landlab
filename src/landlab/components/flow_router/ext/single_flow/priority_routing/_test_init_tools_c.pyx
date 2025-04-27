# distutils: language = c++
# distutils: define_macros=NPY_NO_DEPRECATED_API=NPY_1_7_API_VERSION

import numpy as np

cimport numpy as cnp

from numpy.testing import assert_array_equal

from landlab.components.flow_router.ext.single_flow.priority_routing.init_tools import (
    _get_start_end_indexes_in_sorted_array,
)


def test_get_start_end_indexes_in_sorted_array():
    cdef cnp.int64_t [:] sorted_array = np.int64([0, 0, 0, 1, 2, 2, 3, 3, 3, 3])
    cdef cnp.int64_t n_values = 4, max_value = 9
    cdef cnp.int64_t [:, :] idx2 = np.empty((2, n_values), dtype=np.int64)

    idx2 = _get_start_end_indexes_in_sorted_array(sorted_array, n_values, max_value)
    assert_array_equal(idx2, np.int64([[0, 3, 4, 6], [2, 3, 5, 9]]))

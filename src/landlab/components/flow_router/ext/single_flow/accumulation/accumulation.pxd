# distutils: language = c++

import numpy as np

cimport cython
cimport numpy as cnp
from libcpp cimport bool


cdef cnp.int64_t _add_to_upstream_ordered_nodes(cnp.int64_t receiver_id, cnp.int64_t node_idx_in_stack,
    cnp.int64_t [:] upstream_ordered_nodes, cnp.int64_t [:] donors_start_indexes, cnp.int64_t [:] donors)

cdef void _calc_upstream_order_for_nodes_c(cnp.int64_t[:] base_level_and_closed_nodes,
    cnp.int64_t[:] upstream_ordered_nodes, cnp.int64_t[:] donors_start_indexes, cnp.int64_t[:] donors)

# distutils: language = c++

import numpy as np
cimport numpy as cnp
cimport cython
from libcpp cimport bool

cdef cnp.int_t _add_to_upstream_ordered_nodes(cnp.int_t receiver_id, cnp.int_t node_idx_in_stack, 
    cnp.int_t [:] upstream_ordered_nodes, cnp.int_t [:] donors_start_indexes, cnp.int_t [:] donors)

cdef void _calc_upstream_order_for_nodes_c(cnp.int_t[:] base_level_and_closed_nodes,  
    cnp.int_t[:] upstream_ordered_nodes, cnp.int_t[:] donors_start_indexes, cnp.int_t[:] donors)
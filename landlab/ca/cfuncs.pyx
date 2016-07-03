"""
Created on Thu Jun 30 12:40:39 2016

@author: gtucker
"""

import numpy as np
cimport numpy as np
cimport cython
from landlab import CORE_NODE as _CORE

DTYPE_INT = np.int
ctypedef np.int_t DTYPE_INT_t

DTYPE_INT8 = np.int8
ctypedef np.int8_t DTYPE_INT8_t


# def froggo(a, b):
#     
#     c = a + b

#def update_node_states_cython(a, b):
#
#    cdef int frog
#
#    frog = a + b


def update_node_states_cython(np.ndarray[DTYPE_INT_t, ndim=1] node_state,
                              np.ndarray[DTYPE_INT8_t, ndim=1] status_at_node,
                              DTYPE_INT_t tail_node, 
                              DTYPE_INT_t head_node,
                              DTYPE_INT_t new_link_state,
                              node_pair):

    # Remember the previous state of each node so we can detect whether the
    # state has changed
    old_tail_node_state = node_state[tail_node]
    old_head_node_state = node_state[head_node]

    # Change to the new states
    if status_at_node[tail_node] == _CORE:
        node_state[tail_node] = node_pair[new_link_state][0]
    if status_at_node[head_node] == _CORE:
        node_state[head_node] = node_pair[new_link_state][1]

    return node_state[tail_node] != old_tail_node_state, \
           node_state[head_node] != old_head_node_state

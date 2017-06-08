import numpy as np
cimport numpy as np
cimport cython
from math import exp

DTYPE_FLOAT = np.double
ctypedef np.double_t DTYPE_FLOAT_t

DTYPE_INT = np.int
ctypedef np.int_t DTYPE_INT_t


def calculate_qs_in(np.ndarray[DTYPE_INT_t, ndim=1] stack_flip_ud,
np.ndarray[DTYPE_INT_t, ndim=1] flow_receivers,
    np.ndarray[DTYPE_FLOAT_t, ndim=1] link_lengths,
    np.ndarray[DTYPE_FLOAT_t, ndim=1] q,
    np.ndarray[DTYPE_FLOAT_t, ndim=1] qs,
    np.ndarray[DTYPE_FLOAT_t, ndim=1] qs_in,
    np.ndarray[DTYPE_FLOAT_t, ndim=1] Es,
    np.ndarray[DTYPE_FLOAT_t, ndim=1] Er,
    DTYPE_FLOAT_t v_s,
    DTYPE_FLOAT_t F_f):
    """Calculate qs.
    
    """
    cdef unsigned int n_nodes = stack_flip_ud.size
    cdef unsigned int j
    cdef unsigned int i
    
    #iterate top to bottom through the stack, calculate qs
    for i in range(n_nodes):
        j = stack_flip_ud[i]
        if q[j] == 0:
            qs[j] = 0
        else:
            qs[j] = (((Es[j]) + (1-F_f) * Er[j]) / \
            (v_s / q[j])) * (1.0 - \
            exp(-link_lengths[j] * v_s / q[j])) + \
            (qs_in[j] * exp(-link_lengths[j] * \
            v_s / q[j]))
        qs_in[flow_receivers[j]] += qs[j]

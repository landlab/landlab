import numpy as np

cimport numpy as np


cdef extern from "math.h":
    double exp(double x) nogil

DTYPE_FLOAT = np.double
ctypedef np.double_t DTYPE_FLOAT_t

DTYPE_INT = int
ctypedef np.int_t DTYPE_INT_t


def calculate_qs_in(np.ndarray[DTYPE_INT_t, ndim=1] stack_flip_ud,
                    np.ndarray[DTYPE_INT_t, ndim=1] flow_receivers,
                    np.ndarray[DTYPE_FLOAT_t, ndim=1] cell_area_at_node,
                    np.ndarray[DTYPE_FLOAT_t, ndim=1] q,
                    np.ndarray[DTYPE_FLOAT_t, ndim=1] qs,
                    np.ndarray[DTYPE_FLOAT_t, ndim=1] qs_in,
                    np.ndarray[DTYPE_FLOAT_t, ndim=1] Es,
                    DTYPE_FLOAT_t v_s,
                    DTYPE_FLOAT_t F_f):

    """Calculate and qs and qs_in."""
    # define internal variables
    cdef unsigned int n_nodes = stack_flip_ud.shape[0]
    cdef unsigned int node_id
    cdef unsigned int i

    # iterate top to bottom through the stack, calculate qs and adjust qs_in
    for i in range(n_nodes):

        # choose the node id
        node_id = stack_flip_ud[i]

        # If q at current node is greather than zero, calculate qs based on a
        # local analytical solution. This local analytical solution depends on
        # qs_in, the sediment flux coming into the node from upstream (hence
        # the upstream to downstream node ordering).

        # Because calculation of qs requires qs_in, this operation must be done
        # in an upstream to downstream loop, and cannot be vectorized.
        #
        # there is water flux (q) and this node is not a pit then calculate qs.

        if q[node_id] > 0 and (flow_receivers[node_id] != node_id):
            qs[node_id] = ((qs_in[node_id]
                            + ((1.0 - F_f) * Es[node_id]) * cell_area_at_node[node_id])
                           / (1.0 + (v_s * cell_area_at_node[node_id] / (q[node_id]))))

            # finally, add this node's qs to recieiving nodes qs_in.
            # if qs[node_id] == 0, then there is no need for this line to be
            # evaluated.
            qs_in[flow_receivers[node_id]] += qs[node_id]

        else:
            # if q at the current node is zero, set qs at that node is zero.
            qs[node_id] = 0

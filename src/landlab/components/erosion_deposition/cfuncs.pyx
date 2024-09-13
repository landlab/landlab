cimport cython

# https://cython.readthedocs.io/en/stable/src/userguide/fusedtypes.html
ctypedef fused id_t:
    cython.integral
    long long


cdef extern from "math.h":
    double exp(double x) nogil


@cython.boundscheck(False)
@cython.wraparound(False)
cpdef void calculate_qs_in(
    const id_t [:] stack_flip_ud,
    const id_t [:] flow_receivers,
    const cython.floating [:] cell_area_at_node,
    const cython.floating [:] q,
    cython.floating [:] qs,
    cython.floating [:] qs_in,
    const cython.floating [:] Es,
    const cython.floating [:] v_s,
    const double F_f,
) noexcept nogil:
    """Calculate and qs and qs_in."""
    cdef unsigned int n_nodes = len(stack_flip_ud)
    cdef unsigned int node
    cdef unsigned int i

    # iterate top to bottom through the stack, calculate qs and adjust qs_in
    for i in range(n_nodes):

        # choose the node id
        node = stack_flip_ud[i]

        # If q at current node is greather than zero, calculate qs based on a
        # local analytical solution. This local analytical solution depends on
        # qs_in, the sediment flux coming into the node from upstream (hence
        # the upstream to downstream node ordering).

        # Because calculation of qs requires qs_in, this operation must be done
        # in an upstream to downstream loop, and cannot be vectorized.
        #
        # there is water flux (q) and this node is not a pit then calculate qs.

        if q[node] > 0 and flow_receivers[node] != node:
            qs[node] = (
                (
                    qs_in[node]
                    + (1.0 - F_f) * Es[node] * cell_area_at_node[node]
                ) / (1.0 + v_s[node] * cell_area_at_node[node] / q[node])
            )

            # finally, add this node's qs to recieiving nodes qs_in.
            # if qs[node] == 0, then there is no need for this line to be
            # evaluated.
            qs_in[flow_receivers[node]] += qs[node]

        else:
            # if q at the current node is zero, set qs at that node is zero.
            qs[node] = 0

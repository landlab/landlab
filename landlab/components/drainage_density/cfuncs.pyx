import numpy as np

cimport numpy as np
cimport cython


DTYPE = np.int
ctypedef np.int_t DTYPE_INT_t
DTYPE_FLOAT = np.double
ctypedef np.double_t DTYPE_FLOAT_t


@cython.boundscheck(False)
def _calc_dists_to_channel(np.ndarray[np.uint8_t, ndim=1] ch_network,
                           np.ndarray[DTYPE_INT_t, ndim=1] flow_receivers,
                           np.ndarray[DTYPE_INT_t, ndim=1] upstream_order,
                           np.ndarray[DTYPE_FLOAT_t, ndim=1] link_lengths,
                           np.ndarray[DTYPE_INT_t, ndim=1] stack_links,
                           np.ndarray[DTYPE_FLOAT_t, ndim=1] dist_to_ch,
                           DTYPE_INT_t num_nodes):
    """Calculate distance to nearest channel.

    Calculate the distances to the closest channel node for all nodes in the
    grid.

    Parameters
    ----------
    ch_network : node array
        integer logical map of which nodes contain channels.
    flow_receivers : node array
        ID of the next downstream node.
    link_lengths : num_d8_links-length array of floats
        The length of all links on the grid, including diagonals if present.
    stack_links : node array
        The ID of the link that leads to the downstream node.
    dists_to_ch : number_of_nodes-length array of floats
        The output array; the distance to the nearest channel node.
    num_nodes : int
        The number of nodes.
    """
    cdef DTYPE_INT_t node
    cdef DTYPE_INT_t node_iter
    cdef DTYPE_INT_t flag
    cdef DTYPE_FLOAT_t distance

    for i in range(num_nodes):
        node = upstream_order[i]
        if ch_network[node] == 1:
            # ^we're standing in a channel
            dist_to_ch[node] = 0.
        else:
            distance = 0
            flag = 0
            node_iter = node
            while flag == 0:  # as long as we haven't hit a channel yet...
                if flow_receivers[node_iter] == node_iter:
                    # ^if no flow receiver (boundary probably)
                    ch_network[node_iter] = 1
                    # ^convince the node it's a channel
                else:
                    distance += link_lengths[stack_links[node_iter]]
                    node_iter = flow_receivers[node_iter]
                if ch_network[node_iter] == 1:
                    # ^we've hit a channel
                    dist_to_ch[node] = distance
                    # ^save distance to channel
                    flag = 1

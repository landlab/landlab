cimport cython
from libc.stdint cimport uint8_t

# https://cython.readthedocs.io/en/stable/src/userguide/fusedtypes.html
ctypedef fused id_t:
    cython.integral
    long long


@cython.boundscheck(False)
def _calc_dists_to_channel(
    uint8_t [:] ch_network,
    id_t [:] flow_receivers,
    id_t [:] upstream_order,
    const cython.floating [:] link_lengths,
    id_t [:] stack_links,
    cython.floating [:] dist_to_ch,
    long num_nodes,
):
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
    cdef long node
    cdef long node_iter
    cdef long flag
    cdef double distance

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

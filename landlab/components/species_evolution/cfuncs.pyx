import numpy as np
cimport numpy as np
cimport cython


DTYPE_INT = np.int
ctypedef np.int_t DTYPE_INT_t


def _get_watershed_mask(np.ndarray[DTYPE_INT_t, ndim=1] grid_nodes,
                        np.ndarray[DTYPE_INT_t, ndim=1] flow_receivers,
                        DTYPE_INT_t outlet_id,
                        np.ndarray[DTYPE_INT_t, ndim=1] watershed_mask):
    """

    """
    cdef DTYPE_INT_t node
    cdef DTYPE_INT_t receiver_node

    for node in grid_nodes:
        # Follow flow path of each node.
        receiver_node = flow_receivers[node]
        outlet_not_found = True

        while outlet_not_found:
            node_flows_to_outlet = any([receiver_node == outlet_id,
                                        flow_receivers[receiver_node] ==
                                        outlet_id])
            node_is_outlet = node == outlet_id

            if node_flows_to_outlet or node_is_outlet:
                watershed_mask[node] = True
                outlet_not_found = False

            else:
                receiver_node = flow_receivers[receiver_node]

                if receiver_node == flow_receivers[receiver_node]:
                    # Receiver_node is a pit.
                    outlet_not_found = False

    return watershed_mask

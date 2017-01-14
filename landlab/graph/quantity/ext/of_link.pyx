import numpy as np
cimport numpy as np
cimport cython

DTYPE = np.int
ctypedef np.int_t DTYPE_t


@cython.boundscheck(False)
def calc_midpoint_of_link(np.ndarray[DTYPE_t, ndim=2] nodes_at_link,
                          np.ndarray[np.float_t, ndim=1] x_of_node,
                          np.ndarray[np.float_t, ndim=1] y_of_node,
                          np.ndarray[np.float_t, ndim=2] xy_of_link):
    cdef int link
    cdef int n_links = nodes_at_link.shape[0]

    for link in range(n_links):
        link_tail = nodes_at_link[link][0]
        link_head = nodes_at_link[link][1]

        xy_of_link[link][0] = (x_of_node[link_tail] +
                               x_of_node[link_head]) * .5
        xy_of_link[link][1] = (y_of_node[link_tail] +
                               y_of_node[link_head]) * .5

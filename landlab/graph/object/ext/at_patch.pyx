import numpy as np
cimport numpy as np
cimport cython

from libc.stdlib cimport malloc, free, qsort


from ...ext.argsort cimport unique_int


@cython.boundscheck(False)
def get_nodes_at_patch(np.ndarray[long, ndim=2, mode="c"] links_at_patch,
                       np.ndarray[long, ndim=1, mode="c"] nodes_at_link,
                       np.ndarray[long, ndim=2, mode="c"] nodes_at_patch):
    cdef int n_patches = links_at_patch.shape[0]
    cdef int patch
    cdef int link
    cdef int i
    cdef int n_nodes
    cdef int * unique <int *>malloc(2 * links_at_patch.shape[1] * sizeof(int))
    cdef int * all_nodes <int *>malloc(2 * links_at_patch.shape[1] *
                                       sizeof(int))

    try:
        for patch in range(n_patches):
            i = 0
            n_nodes = 0
            while 1:
                link = links_at_patch[patch, i]
                if link == -1:
                  break

                all_nodes[n_nodes] = nodes_at_link[link][0]
                all_nodes[n_nodes + 1] = nodes_at_link[link][1]

                n_nodes += 2
                i += 1

            n_unique = unique_int(all_nodes, n_nodes, unique)
            for i in range(n_unique):
                nodes_at_patch[patch, i] = unique[i]
      finally:
          free(all_nodes)
          free(unique)

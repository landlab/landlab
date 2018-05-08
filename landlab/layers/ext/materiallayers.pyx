import numpy as np
cimport numpy as np
cimport cython


@cython.boundscheck(False)
def just_erode(np.ndarray[np.float_t, ndim=2] layers,
               int n_layers,
               np.ndarray[np.float_t, ndim=1] dz):
    cdef int n_stacks = layers.shape[1]
    cdef int top_ind = n_layers
    cdef int col
    cdef int layer
    cdef double removed
    cdef double amount_to_remove

    for col in range(n_stacks):

        amount_to_remove = - dz[col]
        removed = 0.
        for layer in range(top_ind - 1, -1, -1):
            removed += layers[layer, col]
            layers[layer, col] = 0.
            if removed > amount_to_remove:
                layers[layer, col] = removed - amount_to_remove
                break

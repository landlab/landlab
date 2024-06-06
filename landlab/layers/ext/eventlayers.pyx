cimport cython

ctypedef fused id_t:
    cython.integral
    long long


@cython.boundscheck(False)
@cython.wraparound(False)
def deposit_or_erode(
    cython.floating [:, :] layers,
    long n_layers,
    const cython.floating [:] dz,
):
    cdef int n_stacks = layers.shape[1]
    cdef int top_ind = n_layers - 1
    cdef int col
    cdef int layer
    cdef double removed
    cdef double amount_to_remove

    with nogil:
        for col in range(n_stacks):
            if dz[col] >= 0.:
                layers[top_ind, col] += dz[col]
            else:
                amount_to_remove = - dz[col]
                removed = 0.
                for layer in range(top_ind, -1, -1):
                    removed += layers[layer, col]
                    layers[layer, col] = 0.
                    if removed > amount_to_remove:
                        layers[layer, col] = removed - amount_to_remove
                        break


@cython.boundscheck(False)
@cython.wraparound(False)
def get_surface_index(
    const cython.floating [:, :] layers,
    long n_layers,
    id_t [:] surface_index,
):
    cdef int n_stacks = layers.shape[1]
    cdef int top_ind = n_layers
    cdef int col
    cdef int layer

    with nogil:
        for col in range(n_stacks):
            for layer in range(top_ind - 1, -1, -1):
                if layers[layer, col] > 0:
                    surface_index[col] = layer
                    break

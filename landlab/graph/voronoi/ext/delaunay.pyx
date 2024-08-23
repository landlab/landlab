cimport cython
from libc.stdlib cimport free
from libc.stdlib cimport malloc

ctypedef fused id_t:
    cython.integral
    long long


@cython.boundscheck(False)
@cython.wraparound(False)
def remove_patches(
    id_t [:, :] links_at_patch,
    id_t [:] patches_to_remove,
):
    cdef int n_patches = links_at_patch.shape[0]
    cdef int max_links = links_at_patch.shape[1]
    cdef int patch
    cdef int n

    new_patch = 0
    for patch in range(n_patches):
        if patch in patches_to_remove:
            pass
        else:
            for n in range(max_links):
                links_at_patch[new_patch, n] = links_at_patch[patch, n]
            patch += 1


@cython.boundscheck(False)
@cython.wraparound(False)
def remove_tris(
    id_t [:, :] nodes_at_tri,
    id_t [:, :] neighbors_at_tri,
    const id_t [:] bad_tris,
):
    cdef int n_tris = nodes_at_tri.shape[0]
    cdef int tri
    cdef int patch
    cdef int *patch_at_tri = <int *>malloc(n_tris * sizeof(int))
    cdef int n
    cdef int old

    try:
        patch = 0
        for tri in range(n_tris):
            if tri in bad_tris:
                patch_at_tri[tri] = -1
            else:
                for n in range(3):
                    nodes_at_tri[patch, n] = nodes_at_tri[tri, n]
                    neighbors_at_tri[patch, n] = neighbors_at_tri[tri, n]
                patch_at_tri[tri] = patch
                patch += 1

        for tri in range(n_tris):
            for n in range(3):
                old = neighbors_at_tri[tri, n]
                if old >= 0:
                    neighbors_at_tri[tri, n] = patch_at_tri[old]
                else:
                    neighbors_at_tri[tri, n] = -1
    finally:
        free(patch_at_tri)


@cython.boundscheck(False)
@cython.wraparound(False)
def _setup_links_at_patch(
    id_t [:, :] nodes_at_patch,
    id_t [:, :] tri_neighbors,
    id_t [:, :] nodes_at_link,
    id_t [:, :] links_at_patch,
):
    cdef int i
    cdef int link
    cdef int neighbor
    cdef int n_patches = len(nodes_at_patch)
    cdef int *tri_done = <int *>malloc(n_patches * sizeof(int))
    cdef int *links_per_patch = <int *>malloc(n_patches * sizeof(int))

    if not tri_done or not links_per_patch:
        raise MemoryError(
            "unable to allocate {bytes} bytes".format(bytes=n_patches * sizeof(int))
        )

    try:
        for tri in range(n_patches):
            tri_done[tri] = 0
            links_per_patch[tri] = 0

        link = 0
        for tri in range(n_patches):
            for i in (0, 1, 2):
                neighbor = tri_neighbors[tri, i]

                if neighbor == -1 or not tri_done[neighbor]:
                    nodes_at_link[link, 0] = nodes_at_patch[tri, (i + 1) % 3]
                    nodes_at_link[link, 1] = nodes_at_patch[tri, (i + 2) % 3]

                    links_at_patch[tri, links_per_patch[tri]] = link
                    links_per_patch[tri] += 1

                    if neighbor >= 0:
                        links_at_patch[neighbor, links_per_patch[neighbor]] = link
                        links_per_patch[neighbor] += 1
                        tri_done[tri] = True

                    link += 1
    finally:
        free(links_per_patch)
        free(tri_done)

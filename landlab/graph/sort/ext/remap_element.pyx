import numpy as np
cimport numpy as np
cimport cython

from libc.stdlib cimport malloc, free


cdef extern from "math.h":
    double atan2(double y, double x) nogil


from .spoke_sort import sort_spokes_at_wheel
from .argsort cimport argsort_int

DTYPE = np.int
ctypedef np.int_t DTYPE_t


@cython.boundscheck(False)
def reverse_one_to_one(np.ndarray[DTYPE_t, ndim=1] mapping,
                       np.ndarray[DTYPE_t, ndim=1] out):
    cdef int n_elements = mapping.size
    cdef int index
    cdef int id_

    for index in range(n_elements):
        id_ = mapping[index]
        if id_ >= 0:
          out[id_] = index


@cython.boundscheck(False)
def reverse_one_to_many(np.ndarray[DTYPE_t, ndim=2] mapping,
                        np.ndarray[DTYPE_t, ndim=2] out):
    cdef int n_elements = mapping.shape[0]
    cdef int n_cols = mapping.shape[1]
    cdef int out_rows = out.shape[0]
    cdef int index
    cdef int id_
    cdef int *count = <int *>malloc(out_rows * sizeof(int))

    try:
        for index in range(out_rows):
            count[index] = 0

        for index in range(n_elements):
            for col in range(n_cols):
                id_ = mapping[index, col]
                if id_ >= 0:
                    out[id_, count[id_]] = index
                    count[id_] += 1
    finally:
        free(count)


@cython.boundscheck(False)
def remap_graph_element(np.ndarray[DTYPE_t, ndim=1] elements,
                        np.ndarray[DTYPE_t, ndim=1] old_to_new):
    """Remap elements in an array in place.

    Parameters
    ----------
    elements : ndarray of int
        Identifiers of elements.
    old_to_new : ndarray of int
        Mapping from the old identifier to the new identifier.
    """
    cdef int n_elements = elements.size
    cdef int i

    for i in range(n_elements):
      elements[i] = old_to_new[elements[i]]


@cython.boundscheck(False)
def remap_graph_element_ignore(np.ndarray[DTYPE_t, ndim=1] elements,
                               np.ndarray[DTYPE_t, ndim=1] old_to_new,
                               DTYPE_t bad_val):
    """Remap elements in an array in place, ignoring bad values.

    Parameters
    ----------
    elements : ndarray of int
        Identifiers of elements.
    old_to_new : ndarray of int
        Mapping from the old identifier to the new identifier.
    bad_val : int
        Ignore values in the input array when remapping.
    """
    cdef int n_elements = elements.size
    cdef int i

    for i in range(n_elements):
        if elements[i] != bad_val:
            elements[i] = old_to_new[elements[i]]


@cython.boundscheck(False)
def reorder_patches(np.ndarray[DTYPE_t, ndim=1] links_at_patch,
                    np.ndarray[DTYPE_t, ndim=1] offset_to_patch,
                    np.ndarray[DTYPE_t, ndim=1] sorted_patches):
    cdef int i
    cdef int patch
    cdef int offset
    cdef int n_links
    cdef int n_patches = len(sorted_patches)
    cdef int *new_offset = <int *>malloc(len(offset_to_patch) * sizeof(int))
    cdef int *new_patches = <int *>malloc(len(links_at_patch) * sizeof(int))

    try:
        new_offset[0] = 0
        for patch in range(n_patches):
            offset = offset_to_patch[sorted_patches[patch]]
            n_links = offset_to_patch[sorted_patches[patch] + 1] - offset

            new_offset[patch + 1] = new_offset[patch] + n_links
            for i in range(n_links):
                new_patches[new_offset[patch] + i] = links_at_patch[offset + i]

        for i in range(len(links_at_patch)):
            links_at_patch[i] = new_patches[i]
        for i in range(len(offset_to_patch)):
            offset_to_patch[i] = new_offset[i]
    finally:
        free(new_offset)
        free(new_patches)


@cython.boundscheck(False)
def calc_center_of_patch(np.ndarray[DTYPE_t, ndim=1] links_at_patch,
                         np.ndarray[DTYPE_t, ndim=1] offset_to_patch,
                         np.ndarray[np.float_t, ndim=2] xy_at_link,
                         np.ndarray[np.float_t, ndim=2] xy_at_patch):
    cdef int patch
    cdef int link
    cdef int i
    cdef int offset
    cdef int n_links
    cdef int n_patches = len(xy_at_patch)
    cdef float x
    cdef float y

    for patch in range(n_patches):
        offset = offset_to_patch[patch]
        n_links = offset_to_patch[patch + 1] - offset
        x = 0.
        y = 0.
        for i in range(offset, offset + n_links):
            link = links_at_patch[i]
            x += xy_at_link[link, 0]
            y += xy_at_link[link, 1]
        xy_at_patch[patch, 0] = x / n_links
        xy_at_patch[patch, 1] = y / n_links


@cython.boundscheck(False)
def reorder_links_at_patch(np.ndarray[DTYPE_t, ndim=1] links_at_patch,
                           np.ndarray[DTYPE_t, ndim=1] offset_to_patch,
                           np.ndarray[np.float_t, ndim=2] xy_of_link):
    cdef int n_patches = len(offset_to_patch) - 1

    xy_of_patch = np.empty((n_patches, 2), dtype=float)
    calc_center_of_patch(links_at_patch, offset_to_patch, xy_of_link,
                         xy_of_patch)
    sort_spokes_at_wheel(links_at_patch, offset_to_patch, xy_of_patch,
                         xy_of_link)


cdef _argsort_links(long * links, int n_links, long * nodes, long * ordered):
    cdef int n_nodes = 2 * n_links
    cdef int * index = <int *>malloc(n_nodes * sizeof(int))

    try:
        argsort_int(nodes, n_nodes, index)

        i = 0
        for link in range(n_links):
            print index[i]
            # if nodes[index[i]] != nodes[index[i + 1]]:
            #     raise ValueError('open patch')

            ordered[i / 2] = index[i] / 2
            i += 2

            # if nodes[index[i]] == nodes[index[i - 1]]:
            #     raise ValueError('triple-point')
    finally:
        free(index)


@cython.boundscheck(False)
def connect_links(np.ndarray[long, ndim=1, mode="c"] links,
                  np.ndarray[long, ndim=2] nodes_at_link):
    cdef long n_links = links.size
    cdef long * nodes = <long *>malloc(2 * n_links * sizeof(long))
    cdef long * ordered = <long *>malloc(n_links * sizeof(long))
    cdef long * buff = <long *>malloc(n_links * sizeof(long))
    cdef long node
    cdef long link

    try:
        node = 0
        for link in range(n_links):
            nodes[node] = nodes_at_link[link, 0]
            nodes[node + 1] = nodes_at_link[link, 1]
            node += 2

        _argsort_links(&links[0], n_links, nodes, ordered)

        for link in range(n_links):
            buff[link] = links[ordered[link]]

        for link in range(n_links):
            links[link] = buff[link]

    finally:
        free(buff)
        free(ordered)
        free(nodes)


cdef reverse_order(long * array, long size):
    cdef long i
    cdef long temp

    for i in range(size / 2):
        temp = array[i]
        array[i] = array[(size - 1) - i]
        array[(size - 1) - i] = temp

        # array[i], array[size - 1] = array[size - 1], array[i]


@cython.boundscheck(False)
def reverse_element_order(np.ndarray[long, ndim=2] links_at_patch,
                          np.ndarray[long, ndim=1] patches):
    cdef long n_patches = patches.shape[0]
    cdef long max_links = links_at_patch.shape[1]
    cdef long patch
    cdef long n
    cdef long i

    for i in range(n_patches):
        patch = patches[i]
        # for n in range(1, max_links):
        #     if links_at_patch[patch, n] == -1:
        #         break
        # reverse_order(&links_at_patch[patch, 1], n - 1)

        n = 1
        while n < max_links:
            if links_at_patch[patch, n] == -1:
                break
            n += 1
        reverse_order(&links_at_patch[patch, 1], n - 1)

@cython.boundscheck(False)
def get_angle_of_link(np.ndarray[DTYPE_t, ndim=2] nodes_at_link,
                      np.ndarray[np.float_t, ndim=2] xy_of_node,
                      np.ndarray[np.float_t, ndim=1] angle_of_link):
    cdef int link
    cdef float link_tail_x
    cdef float link_tail_y
    cdef float link_head_x
    cdef float link_head_y
    cdef int n_links = nodes_at_link.shape[0]

    for link in range(n_links):
        link_tail_x = xy_of_node[nodes_at_link[link][0]][0]
        link_tail_y = xy_of_node[nodes_at_link[link][0]][1]
        link_head_x = xy_of_node[nodes_at_link[link][1]][0]
        link_head_y = xy_of_node[nodes_at_link[link][1]][1]

        angle_of_link[link] = atan2(link_head_y - link_tail_y,
                                    link_head_x - link_head_y)


@cython.boundscheck(False)
def reorient_links(np.ndarray[DTYPE_t, ndim=2] nodes_at_link,
                   np.ndarray[DTYPE_t, ndim=1] xy_of_node):
    """Reorient links to point up and to the right.

    Parameters
    ----------
    nodes_at_link : ndarray of int, shape `(n_nodes, 2)`
        Identifier for node at link tail and head.
    xy_of_node : ndarray of float, shape `(n_nodes, 2)`
        Coordinate of node as `(x, y)`.
    """
    cdef int link
    cdef int temp
    cdef double angle
    cdef int n_links = nodes_at_link.shape[0]
    cdef double minus_45 = - np.pi * .25
    cdef double plus_135 = np.pi * .75

    angle_of_link = np.empty(n_links, dtype=n_links)
    get_angle_of_link(nodes_at_link, xy_of_node, angle_of_link)

    for link in range(n_links):
        angle = angle_of_link[link]
        if angle < minus_45 or angle > plus_135:
            temp = nodes_at_link[link, 0]
            nodes_at_link[link, 0] = nodes_at_link[link, 1]
            nodes_at_link[link, 1] = temp

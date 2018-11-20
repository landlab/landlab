import numpy as np
cimport numpy as np
cimport cython

from libc.stdlib cimport malloc, free

DTYPE = np.int
ctypedef np.int_t DTYPE_t


@cython.boundscheck(False)
def calc_area_at_patch(np.ndarray[long, ndim=2] nodes_at_patch,
                       np.ndarray[np.float_t, ndim=1] x_of_node,
                       np.ndarray[np.float_t, ndim=1] y_of_node,
                       np.ndarray[np.float_t, ndim=1] out):
    cdef long n_patches = nodes_at_patch.shape[0]
    cdef long n_vertices = nodes_at_patch.shape[1]
    cdef long n

    for n in range(n_patches):
      out[n] = calc_area_of_patch(&nodes_at_patch[n, 0], n_vertices,
                                  &x_of_node[0], &y_of_node[0])


cdef calc_area_of_patch(long * nodes_at_patch, long n_vertices,
                        double * x_of_node, double * y_of_node):
    cdef int n
    cdef int node
    cdef double * x_of_vertex = <double *>malloc(n_vertices * sizeof(double))
    cdef double * y_of_vertex = <double *>malloc(n_vertices * sizeof(double))

    try:
        for n in range(n_vertices):
            node = nodes_at_patch[n]
            if node == -1:
                n -= 1
                break
            x_of_vertex[n] = x_of_node[node]
            y_of_vertex[n] = y_of_node[node]

        return calc_area_of_polygon(x_of_vertex, y_of_vertex, n + 1)
    finally:
        free(y_of_vertex)
        free(x_of_vertex)


@cython.boundscheck(False)
def calc_centroid_at_patch(np.ndarray[long, ndim=2] nodes_at_patch,
                           np.ndarray[np.float_t, ndim=1] x_of_node,
                           np.ndarray[np.float_t, ndim=1] y_of_node,
                           np.ndarray[np.float_t, ndim=2] out):
    cdef long n_patches = nodes_at_patch.shape[0]
    cdef long n_vertices = nodes_at_patch.shape[1]
    cdef long n

    for n in range(n_patches):
        calc_centroid_of_patch(&nodes_at_patch[n, 0], n_vertices,
                               &x_of_node[0], &y_of_node[0],
                               &out[n, 0])


cdef calc_centroid_of_patch(long * nodes_at_patch, long n_vertices,
                            double * x_of_node, double * y_of_node, double * out):
    cdef int n
    cdef int node
    cdef double * x = <double *>malloc(n_vertices * sizeof(double))
    cdef double * y = <double *>malloc(n_vertices * sizeof(double))

    try:
        for n in range(n_vertices):
            node = nodes_at_patch[n]
            if node == -1:
                n -= 1
                break
            x[n] = x_of_node[node]
            y[n] = y_of_node[node]
        calc_centroid_of_polygon(x, y, n + 1, out)
    finally:
        free(y)
        free(x)


cdef calc_centroid_of_polygon(double * x, double * y, long n_vertices,
                              double * out):
    cdef double x_of_centroid = 0.
    cdef double y_of_centroid = 0.
    cdef double area = calc_area_of_polygon(x, y, n_vertices)
    cdef double c
    cdef int n

    c = x[n_vertices - 1] * y[0] - x[0] * y[n_vertices - 1]
    x_of_centroid = (x[n_vertices - 1] + x[0]) * c
    y_of_centroid = (y[n_vertices - 1] + y[0]) * c

    for n in range(n_vertices - 1):
        c = x[n] * y[n + 1] - x[n + 1] * y[n]
        x_of_centroid += (x[n] + x[n + 1]) * c
        y_of_centroid += (y[n] + y[n + 1]) * c

    # c = x[n_vertices - 1] * y[0] - x[0] - y[n_vertices - 1]
    # x_of_centroid += (x[n_vertices - 1] + x[0]) * c
    # y_of_centroid += (y[n_vertices - 1] + y[0]) * c

    x_of_centroid /= 6. * area
    y_of_centroid /= 6. * area

    out[0] = x_of_centroid
    out[1] = y_of_centroid


cdef calc_area_of_polygon(double * x, double * y, long n_vertices):
    cdef double area = 0.
    cdef int n

    for n in range(n_vertices - 1):
        area += x[n] * y[n + 1]
        area -= x[n + 1] * y[n]
    area += x[n_vertices - 1] * y[0]
    area -= x[0] * y[n_vertices - 1]
    area /= 2.

    return area

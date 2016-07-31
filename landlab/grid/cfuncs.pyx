import numpy as np
cimport numpy as np
cimport cython


DTYPE = np.int
ctypedef np.int_t DTYPE_INT_t
DTYPE_FLOAT = np.double
ctypedef np.double_t DTYPE_FLOAT_t


@cython.boundscheck(False)
def find_rows_containing_ID(np.ndarray[DTYPE_INT_t, ndim=2] input_array,
                            np.ndarray[DTYPE_INT_t, ndim=2] out):
    """
    Record the row in which ID appears in in_array, indexed by ID.

    Parameters
    ----------
    input_array : 2d array of ints
        The input array. The function will report in which rows it finds
        each ID.
    out : 2d array of ints, num_IDs x num_rows_in_input_array
        The output array. The number of rows sets the number of IDs to search
        for in input_array.
    """

    cdef int nrows = out.shape[0]
    cdef int i
    cdef np.ndarray[DTYPE_INT_t, ndim=2] contains_ID = np.empty_like(
        input_array)

    for i in range(nrows):
        contains_ID = np.equal(input_array, i).astype(int)
        out[i, :] = contains_ID.sum(axis=1)


# @cython.boundscheck(False)
# def _create    


@cython.boundscheck(False)
def _argsort_points_by_x_then_y(np.ndarray[DTYPE_INT_t, ndim=1] pts,
                                np.ndarray[DTYPE_INT_t, ndim=1] out):

    cdef np.ndarray[DTYPE_INT_t, ndim=1] a
    cdef np.ndarray[DTYPE_INT_t, ndim=1] b

    a = pts[:, 0].argsort(kind='mergesort')
    b = pts[a, 1].argsort(kind='mergesort')
    out[:] = a[b]


@cython.boundscheck(False)
def _anticlockwise_argsort_points(np.ndarray[DTYPE_INT_t, ndim=2] pts,
                                  np.ndarray[DTYPE_INT_t, ndim=1] out):

    cdef np.ndarray[DTYPE_FLOAT_t, ndim=1] midpt = pts.mean(axis=0)
    cdef np.ndarray[DTYPE_FLOAT_t, ndim=1] theta = np.empty(
        pts.shape[0], dtype=float)
    cdef double twopi = 2.*np.pi

    theta[:] = np.arctan2(pts[:, 1] - midpt[1], pts[:, 0] - midpt[0])
    theta[:] = theta % twopi
    out[:] = np.argsort(theta)


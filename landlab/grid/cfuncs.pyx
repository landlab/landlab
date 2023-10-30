import numpy as np

cimport cython
cimport numpy as np

DTYPE = int
ctypedef np.int_t DTYPE_INT_t
DTYPE_INTP = np.intp
ctypedef np.intp_t DTYPE_INTP_t
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


@cython.boundscheck(False)
def create_patches_at_element(
        np.ndarray[DTYPE_INT_t, ndim=2] elements_at_patch,
        int number_of_elements, np.ndarray[DTYPE_INT_t, ndim=2] out):
    """
    """
    cdef int i
    cdef np.ndarray[DTYPE_INT_t, ndim=2] element_with_value = np.empty_like(
        (elements_at_patch), dtype=int)
    cdef np.ndarray[DTYPE_INTP_t, ndim=1] patches_with_element
    cdef int num_elements_here

    for i in range(number_of_elements):
        np.equal(elements_at_patch, i, out=element_with_value)
        patches_with_element = np.argwhere(element_with_value)[:, 0]
        num_elements_here = patches_with_element.shape[0]
        out[i, :num_elements_here] = patches_with_element


@cython.boundscheck(False)
def create_links_at_patch(np.ndarray[DTYPE_INT_t, ndim=2] nodes_at_patch,
                          np.ndarray[DTYPE_INT_t, ndim=2] links_at_node,
                          int number_of_patches,
                          np.ndarray[DTYPE_INT_t, ndim=2] out):

    cdef int i
    cdef np.ndarray[DTYPE_INT_t, ndim=1] nodes_on_patch = np.empty(
        nodes_at_patch.shape[1], dtype=int)
    cdef np.ndarray[DTYPE_INT_t, ndim=2] links_at_patch_nodes = np.empty(
        (nodes_at_patch.shape[1], links_at_node.shape[1]), dtype=int)
    cdef np.ndarray[DTYPE_INT_t, ndim=1] vals
    cdef np.ndarray[DTYPE_INTP_t, ndim=1] counts
    cdef np.ndarray[DTYPE_INT_t, ndim=1] duplicated_vals

    for i in range(number_of_patches):
        nodes_on_patch = nodes_at_patch[i, :]
        links_at_patch_nodes = links_at_node[nodes_on_patch, :]
        vals, counts = np.unique(links_at_patch_nodes,
                                 return_counts=True)
        duplicated_vals = vals[counts == 2]
        # if len==4, contains a -1; strip it
        out[i, :] = duplicated_vals[-3:]  # out of order, so far


@cython.boundscheck(False)
def _argsort_points_by_x_then_y(np.ndarray[DTYPE_INT_t, ndim=1] pts,
                                np.ndarray[DTYPE_INT_t, ndim=1] out):

    cdef np.ndarray[DTYPE_INT_t, ndim=1] a
    cdef np.ndarray[DTYPE_INT_t, ndim=1] b

    a = pts[:, 0].argsort(kind="mergesort")
    b = pts[a, 1].argsort(kind="mergesort")
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

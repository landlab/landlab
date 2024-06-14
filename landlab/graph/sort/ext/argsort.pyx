import numpy as np

cimport cython
cimport numpy as np

from cython.parallel import prange

from cython cimport view
from libc.stdlib cimport free
from libc.stdlib cimport malloc

ctypedef np.int_t INT_t
ctypedef np.float_t FLOAT_t

cdef struct Sorter:
    INT_t index
    FLOAT_t value


cdef struct IntSorter:
    INT_t index
    INT_t value


cdef extern from "stdlib.h":
    ctypedef void const_void "const void"
    int qsort(void *, size_t, size_t, int(*)(const_void *, const_void *)) nogil
    int	mergesort(void *, size_t, size_t, int (*)(const_void *, const_void *)) nogil


cdef int _compare(const_void *a, const_void *b) noexcept:
    cdef double v = ((<Sorter*>a)).value - ((<Sorter*>b)).value
    if v < 0:
        return -1
    elif v > 0:
        return 1
    else:
        return 0


cdef int _compare_int(const_void *a, const_void *b) noexcept:
    cdef int v = ((<IntSorter*>a)).value - ((<IntSorter*>b)).value
    if v < 0:
        return -1
    elif v > 0:
        return 1
    else:
        return 0


cdef void argsort_flt(
    cython.floating * data, int n_elements, id_t * out
) noexcept nogil:
    cdef Sorter *sorted_struct = <Sorter*>malloc(n_elements * sizeof(Sorter))
    cdef int i

    try:
        # _argsort(data, n_elements, sorted_struct)
        for i in range(n_elements):
            sorted_struct[i].index = i
            sorted_struct[i].value = data[i]

        qsort(<void*> sorted_struct, n_elements, sizeof(Sorter), _compare_int)

        for i in range(n_elements):
            out[i] = sorted_struct[i].index
    finally:
        free(sorted_struct)


cdef void argsort_int(long * data, int n_elements, int * out) noexcept nogil:
    cdef IntSorter *sorted_struct = <IntSorter*>malloc(n_elements * sizeof(IntSorter))
    cdef int i

    try:
        for i in range(n_elements):
            sorted_struct[i].index = i
            sorted_struct[i].value = data[i]

        qsort(<void*> sorted_struct, n_elements, sizeof(IntSorter), _compare_int)

        for i in range(n_elements):
            out[i] = sorted_struct[i].index
    finally:
        free(sorted_struct)


cdef int unique_int(long * data, int n_elements, int * out):
    cdef int * index = <int *>malloc(n_elements * sizeof(int))
    cdef int n_unique

    try:
        argsort_int(data, n_elements, index)

        if data[index[0]] != data[index[1]]:
            out[0] = data[index[0]]
        else:
            out[0] = data[index[1]]

        n_unique = 1
        for i in range(1, n_elements):
            if data[index[i]] != data[index[i - 1]]:
                out[n_unique] = data[index[i]]
                n_unique += 1
    finally:
        free(index)

    return n_unique


@cython.boundscheck(False)
@cython.wraparound(False)
cpdef sort_children_at_parent(
    id_t [:, :] children_at_parent,
    cython.floating [:] value_at_child,
    id_t [:, :] out,
):
    """Sort the children of parents based on child values.

    Parameters
    ----------
    children_at_parent : ndarray of shape (n_parents, max_family_size)
        Matrix where each row is the IDs of the children (i.e. an index
        into the value_at_child array) belonging to a parent. A value
        of -1 indicates the lack of a child. Note that children can
        belong to multiple parents.
    value_at_child : ndarray of shape (n_children,)
        A value for each child to sort on.
    out : ndarray of shape (n_parents, max_family_size)
        Output array that contains the sorted children for each parent.
        Missing children (i.e. -1 values) are pushed to the each of
        each row.
    """
    cdef int n_parents = children_at_parent.shape[0]
    cdef int n_cols = children_at_parent.shape[1]
    cdef int n_children
    cdef int col
    cdef int parent, child
    cdef double [:, :] values = view.array(
        shape=(n_parents, n_cols),
        itemsize=sizeof(double),
        format="d",
        allocate_buffer=True,
    )
    cdef int [:, :] sorted_indices = view.array(
        shape=(n_parents, n_cols),
        itemsize=sizeof(int),
        format="i",
        allocate_buffer=True,
    )
    cdef int [:, :] indices = view.array(
        shape=(n_parents, n_cols),
        itemsize=sizeof(int),
        format="i",
        allocate_buffer=True,
    )

    for parent in prange(n_parents, nogil=True, schedule="static"):
        n_children = 0
        for col in range(n_cols):
            child = children_at_parent[parent, col]
            if child != -1:
                indices[parent, n_children] = child
                values[parent, n_children] = value_at_child[child]
                n_children = n_children + 1

        if n_children > 0:
            argsort_flt(&values[parent, 0], n_children, &sorted_indices[parent, 0])
            for col in range(n_children):
                out[parent, col] = indices[parent, sorted_indices[parent, col]]
            for col in range(n_children, n_cols):
                out[parent, col] = -1


@cython.boundscheck(False)
@cython.wraparound(False)
cpdef sort_id_array(
    const id_t [:, :] id_array,
    const cython.floating [:, :] data,
    id_t [:, :] out,
):
    """Sort rows of an id-array matrix.

    Parameters
    ----------
    children_at_parent : ndarray of shape (n_parents, max_family_size)
        Matrix where each row is the IDs of the children (i.e. an index
        into the value_at_child array) belonging to a parent. A value
        of -1 indicates the lack of a child.
    value_at_child : ndarray of shape (n_children,)
        A value for each child to sort on.
    out : ndarray of shape (n_parents, max_family_size)
        Output array that contains the sorted children for each parent.
        Missing children (i.e. -1 values) are pushed to the each of
        each row.
    """
    cdef int n_rows = data.shape[0]
    cdef int n_cols = data.shape[1]
    cdef int row, col
    cdef int [:, :] sorted_indices = view.array(
        shape=(n_rows, n_cols),
        itemsize=sizeof(int),
        format="i",
        allocate_buffer=True,
    )
    cdef id_t * temp

    for row in prange(n_rows, nogil=True, schedule="static"):
        temp = <id_t*>malloc(sizeof(id_t) * n_cols)
        try:
            for col in range(n_cols):
                temp[col] = id_array[row, col]

            argsort_id_array(data[row, :], id_array[row, :], sorted_indices[row, :])
            for col in range(n_cols):
                if sorted_indices[row, col] != -1:
                    out[row, col] = temp[sorted_indices[row, col]]
                else:
                    out[row, col] = -1
        finally:
            free(temp)


@cython.boundscheck(False)
@cython.wraparound(False)
cpdef void argsort_id_array(
    const cython.floating [:] data,
    const id_t [:] id_array,
    cython.integral [:] out,
) noexcept nogil:
    cdef int n_elements = len(data)
    cdef int i
    cdef int count = 0
    cdef Sorter *order = <Sorter*>malloc(n_elements * sizeof(Sorter))

    try:
        for i in range(n_elements):
            if id_array[i] != -1:
                order[count].index = i
                order[count].value = data[i]
                count = count + 1

        qsort(<void*> order, count, sizeof(Sorter), _compare)

        for i in range(count):
            out[i] = order[i].index
        for i in range(count, n_elements):
            out[i] = -1
    finally:
        free(order)

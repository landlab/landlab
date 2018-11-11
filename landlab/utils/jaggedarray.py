"""The JaggedArray class to store arrays of variable-length arrays.

Examples
--------

Create a JaggedArray that stores link IDs for the links attached to the
nodes of a 3x3 grid.

>>> from landlab.utils.jaggedarray import JaggedArray
>>> links_at_node = JaggedArray([
...     [0, 6],
...     [1, 7, 0],
...     [8, 1],
...     [2, 9, 6],
...     [3, 10, 2, 7],
...     [11, 3, 8],
...     [4, 7],
...     [5, 10, 4],
...     [5, 11]])

Make up some data that provides values at each of the links.

>>> value_at_link = np.arange(12, dtype=float)

Create another JaggedArray. Here we store the values at each of the links
attached to nodes of the grid.

>>> values_at_node = JaggedArray.empty_like(links_at_node, dtype=float)
>>> values_at_node.array[:] = value_at_link[links_at_node.array]

Now operate on the link values for each node.

>>> values_at_node.foreach_row(sum)
array([  6.,   8.,   9.,  17.,  22.,  22.,  11.,  19.,  16.])
>>> values_at_node.foreach_row(min)
array([ 0.,  0.,  1.,  2.,  2.,  3.,  4.,  4.,  5.])
>>> values_at_node.foreach_row(np.ptp)
array([ 6.,  7.,  7.,  7.,  8.,  8.,  3.,  6.,  6.])
"""
import numpy as np
from six.moves import range


def flatten_jagged_array(jagged, dtype=None):
    """Flatten a list of lists.

    Parameters
    ----------
    jagged : array_like of array_like
        An array of arrays of unequal length.

    Returns
    -------
    (data, offset) : (ndarray, ndarray of int)
        A tuple the data, as a flat numpy array, and offsets into that array
        for every item of the original list.

    Examples
    --------
    >>> from landlab.utils.jaggedarray import flatten_jagged_array
    >>> data, offset = flatten_jagged_array([[1, 2], [], [3, 4, 5]], dtype=int)
    >>> data
    array([1, 2, 3, 4, 5])
    >>> offset
    array([0, 2, 2, 5])
    """
    data = np.concatenate(jagged).astype(dtype=dtype)
    # if len(jagged) > 1:
    #     data = np.concatenate(jagged).astype(dtype=dtype)
    # else:
    #     data = np.array(jagged[0]).astype(dtype=dtype)
    items_per_block = np.array([len(block) for block in jagged], dtype=int)

    offset = np.empty(len(items_per_block) + 1, dtype=int)
    offset[0] = 0
    offset[1:] = np.cumsum(items_per_block)

    return data, offset


def unravel(data, offset, out=None, pad=None):
    """Unravel a jagged array.

    Parameters
    ----------
    data : ndarray
        Flattened-array of the data.
    offset : ndarray of int
        Offsets to the start of rows of the jagged array.
    out : ndarray
        Buffer into which to place the unravelled data.
    pad : number
        Value to use to pad rows of the jagged array.

    Returns
    -------
    ndarray
        Matrix that holds the unravelled jagged array.
    """
    from .ext.jaggedarray import unravel

    n_cols = np.diff(offset).max()
    if out is None:
        if pad is None:
            out = np.empty((len(offset) - 1, n_cols), dtype=data.dtype)
        else:
            out = np.full((len(offset) - 1, n_cols), pad, dtype=data.dtype)
    else:
        if pad is not None:
            out.fill(pad)

    unravel(data, offset, out)

    return out


class JaggedArray(object):

    """
    A container for an array of variable-length arrays.

    JaggedArray([row0, row1, ...])
    JaggedArray(values, values_per_row)

    Examples
    --------
    Create a JaggedArray with an array of arrays.

    >>> from landlab.utils.jaggedarray import JaggedArray
    >>> x = JaggedArray([[0, 1, 2], [3, 4]])
    >>> x.array
    array([0, 1, 2, 3, 4])

    Create a JaggedArray as a 1D array and a list or row lengths.

    >>> x = JaggedArray([0, 1, 2, 3, 4], (3, 2))
    >>> x.array
    array([0, 1, 2, 3, 4])
    """

    def __init__(self, *args):
        """
        JaggedArray([row0, row1, ...])
        JaggedArray(values, values_per_row)

        Examples
        --------
        Create a JaggedArray with an array of arrays.

        >>> from landlab.utils.jaggedarray import JaggedArray
        >>> x = JaggedArray([[0, 1, 2], [3, 4]])
        >>> x.array
        array([0, 1, 2, 3, 4])

        >>> x = JaggedArray([[0, 1, 2]])
        >>> x.array
        array([0, 1, 2])
        >>> x.offset
        array([0, 3])

        Create a JaggedArray as a 1D array and a list or row lengths.

        >>> x = JaggedArray([0, 1, 2, 3, 4], (3, 2))
        >>> x.array
        array([0, 1, 2, 3, 4])
        """
        if len(args) == 1:
            values, values_per_row = (
                np.concatenate(args[0]),
                [len(row) for row in args[0]],
            )
            # if len(args[0]) > 1:
            #     values, values_per_row = (np.concatenate(args[0]),
            #                               [len(row) for row in args[0]])
            # else:
            #     values, values_per_row = (np.array(args[0]), [len(args[0])])
        else:
            values, values_per_row = (np.array(args[0]), args[1])

        self._values = values
        self._number_of_rows = len(values_per_row)
        self._offsets = JaggedArray._offsets_from_values_per_row(values_per_row)
        self._offsets.flags["WRITEABLE"] = False

    @property
    def array(self):
        """The jagged array as a 1D array.

        Returns
        -------
        array :
            A view of the underlying 1D array.

        Examples
        --------
        >>> from landlab.utils.jaggedarray import JaggedArray
        >>> x = JaggedArray([[0, 1, 2], [3, 4]])
        >>> x.array
        array([0, 1, 2, 3, 4])

        >>> x.array[0] = 1
        >>> x.array
        array([1, 1, 2, 3, 4])
        """
        return self._values

    @property
    def offset(self):
        """The offsets to rows of a 1D array.

        Returns
        -------
        array :
            Offsets into the underlying 1D array.

        Examples
        --------
        >>> from landlab.utils.jaggedarray import JaggedArray
        >>> x = JaggedArray([[0, 1, 2], [3, 4]])
        >>> x.offset
        array([0, 3, 5])

        From the offsets you can get values for rows of the jagged array.

        >>> x.array[x.offset[0]:x.offset[1]]
        array([0, 1, 2])

        Once the array is created, you can't change the offsets.

        >>> x.offset[0] = 1
        Traceback (most recent call last):
        ValueError: assignment destination is read-only
        """
        return self._offsets

    @property
    def size(self):
        """Number of array elements.

        Returns
        -------
        int :
            Number of values in the array.

        Examples
        --------
        >>> from landlab.utils.jaggedarray import JaggedArray
        >>> x = JaggedArray([[0, 1, 2], [3, 4]])
        >>> x.size
        5
        """
        return len(self._values)

    @property
    def number_of_rows(self):
        """Number of array rows.

        Returns
        -------
        int :
            Number of rows in the array.

        Examples
        --------
        >>> from landlab.utils.jaggedarray import JaggedArray
        >>> x = JaggedArray([[0, 1, 2], [3, 4]])
        >>> x.number_of_rows
        2
        """
        return self._number_of_rows

    @staticmethod
    def _offsets_from_values_per_row(values_per_row):
        """Get offsets into the base array from array lengths.

        Parameters
        ----------
        values_per_row : array of int
            The number of values in each row of the JaggedArray.

        Returns
        -------
        ndarray
            An array of offsets.
        """
        offset = np.empty(len(values_per_row) + 1, dtype=int)
        np.cumsum(values_per_row, out=offset[1:])
        offset[0] = 0
        return offset

    @staticmethod
    def empty_like(jagged, dtype=None):
        """Create a new JaggedArray that is like another one.

        Parameters
        ----------
        jagged : JaggedArray
            A JaggedArray to copy.
        dtype : np.dtype
            The data type of the new JaggedArray.

        Returns
        -------
        JaggedArray
            A new JaggedArray.
        """
        return JaggedArray(
            np.empty_like(jagged.array, dtype=dtype), np.diff(jagged.offset)
        )

    def length_of_row(self, row):
        """Number of values in a given row.

        Parameters
        ----------
        row : int
            Index to a row.

        Returns
        -------
        int :
            Number of values in the row.

        Examples
        --------
        >>> from landlab.utils.jaggedarray import JaggedArray
        >>> x = JaggedArray([[0, 1, 2], [3, 4]])
        >>> x.length_of_row(0)
        3
        >>> x.length_of_row(1)
        2
        """
        return self._offsets[row + 1] - self._offsets[row]

    def row(self, row):
        """Get the values of a row as an array.

        Parameters
        ----------
        row : int
            Index to a row.

        Returns
        -------
        array :
            Values in the row as a slice of the underlying array.

        Examples
        --------
        >>> from landlab.utils.jaggedarray import JaggedArray
        >>> x = JaggedArray([[0, 1, 2], [3, 4]])
        >>> x.row(0)
        array([0, 1, 2])
        >>> x.row(1)
        array([3, 4])

        >>> y = x.row(0)
        >>> y[0] = 1
        >>> x.row(0)
        array([1, 1, 2])
        """
        return self._values[self._offsets[row] : self._offsets[row + 1]]

    def __iter__(self):
        """Iterate over the rows of the array.

        Examples
        --------
        >>> from landlab.utils.jaggedarray import JaggedArray
        >>> x = JaggedArray([[0, 1, 2], [3, 4]])
        >>> for row in x: row
        array([0, 1, 2])
        array([3, 4])
        """
        for row_number in range(self._number_of_rows):
            yield self.row(row_number)

    def foreach_row(self, func, out=None):
        """Apply an operator row-by-row.

        Examples
        --------
        >>> from landlab.utils.jaggedarray import JaggedArray
        >>> x = JaggedArray([[0, 1, 2], [3, 4]])
        >>> x.foreach_row(sum)
        array([3, 7])

        >>> out = np.empty(2, dtype=int)
        >>> x.foreach_row(sum, out=out) is out
        True
        >>> out
        array([3, 7])
        """
        if out is None:
            out = np.empty(self.number_of_rows, dtype=self._values.dtype)

        for (row_number, row) in enumerate(self):
            out[row_number] = func(row)

        return out

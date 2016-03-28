"""Store arrays of variable-length arrays implemented with masked arrays.

Implements a JaggedArray class using numpy masked arrays.

Examples
--------

Create a JaggedArray that stores link IDs for the links attached to the
nodes of a 3x3 grid.

>>> from landlab.utils.jaggedarray_ma import JaggedArray
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
>>> values_at_node.array = value_at_link[links_at_node.array]

Now operate on the link values for each node.

>>> values_at_node.foreach_row(np.sum)
array([  6.,   8.,   9.,  17.,  22.,  22.,  11.,  19.,  16.])
>>> values_at_node.foreach_row(np.min)
array([ 0.,  0.,  1.,  2.,  2.,  3.,  4.,  4.,  5.])
>>> values_at_node.foreach_row(np.ptp)
array([ 6.,  7.,  7.,  7.,  8.,  8.,  3.,  6.,  6.])
"""
import numpy as np
from six.moves import range


class JaggedArray(object):

    """
    A container for an array of variable-length arrays.

    JaggedArray([row0, row1, ...])
    JaggedArray(values, values_per_row)

    Examples
    --------
    Create a JaggedArray with an array of arrays.

    >>> from landlab.utils.jaggedarray_ma import JaggedArray
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

        >>> from landlab.utils.jaggedarray_ma import JaggedArray
        >>> x = JaggedArray([[0, 1, 2], [3, 4]])
        >>> x.array
        array([0, 1, 2, 3, 4])

        Create a JaggedArray as a 1D array and a list or row lengths.

        >>> x = JaggedArray([0, 1, 2, 3, 4], (3, 2))
        >>> x.array
        array([0, 1, 2, 3, 4])
        """
        if len(args) == 1:
            if isinstance(args[0], np.ma.core.MaskedArray):
                mat = args[0]
            else:
                mat = JaggedArray.ma_from_list_of_lists(args[0])
        else:
            mat = JaggedArray.ma_from_flat_array(args[0], args[1])

        self._values = mat
        self._number_of_rows = mat.shape[0]

    @staticmethod
    def ma_from_list_of_lists(rows, dtype=None):
        """Create a masked array from a list of lists.

        Parameters
        ----------
        rows : array_like or array_like
            Rows of the jagged array.
        dtype : np.dtype, optional
            The data type of the new masked array.

        Returns
        -------
        np.masked_array
            A new masked array.
        """
        values_per_row = [len(row) for row in rows]
        mat = np.ma.masked_all((len(rows), max(values_per_row)),
                               dtype=dtype or int)
        for (row_number, row) in enumerate(rows):
            mat[row_number, :len(row)] = row

        return mat

    @staticmethod
    def ma_from_flat_array(array, values_per_row):
        """Create a masked array from a flat array.

        Parameters
        ----------
        array : array_like
            Values of the jagged array.
        values_per_row : array_like of int
            Number of values in each row of the jagged array.

        Returns
        -------
        np.masked_array
            A new masked array.
        """
        array = np.array(array)
        mat = np.ma.masked_all((len(values_per_row), max(values_per_row)),
                               dtype=array.dtype)
        offset = 0
        for row_number in range(mat.shape[0]):
            n_valid = values_per_row[row_number]
            mat[row_number, :n_valid] = array[offset:offset + n_valid]
            offset += n_valid

        return mat

    @property
    def array(self):
        """The jagged array as a 1D array.

        Returns
        -------
        array :
            A view of the underlying 1D array.

        Examples
        --------
        >>> from landlab.utils.jaggedarray_ma import JaggedArray
        >>> x = JaggedArray([[0, 1, 2], [3, 4]])
        >>> x.array
        array([0, 1, 2, 3, 4])

        >>> x.array = np.array([1, 1, 2, 3, 4])
        >>> x.array
        array([1, 1, 2, 3, 4])
        """
        return self._values.compressed()

    @property
    def masked_array(self):
        """The jagged array as a masked array.

        Returns
        -------
        np.masked_array :
            The underlying masked array.
        """
        return self._values

    @array.setter
    def array(self, array):
        """Set the data of the jagged array from a 1D array.

        Parameters
        ----------
        array : array_like
            The new values of the array.
        """
        self._values[~ self._values.mask] = array

    @property
    def size(self):
        """Number of array elements.

        Returns
        -------
        int :
            Number of values in the array.

        Examples
        --------
        >>> from landlab.utils.jaggedarray_ma import JaggedArray
        >>> x = JaggedArray([[0, 1, 2], [3, 4]])
        >>> x.size
        5
        """
        return self.array.size

    @property
    def number_of_rows(self):
        """Number of array rows.

        Returns
        -------
        int :
            Number of rows in the array.

        Examples
        --------
        >>> from landlab.utils.jaggedarray_ma import JaggedArray
        >>> x = JaggedArray([[0, 1, 2], [3, 4]])
        >>> x.number_of_rows == 2
        True
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
        return JaggedArray(np.ma.empty_like(jagged.masked_array, dtype=dtype))

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
        >>> from landlab.utils.jaggedarray_ma import JaggedArray
        >>> x = JaggedArray([[0, 1, 2], [3, 4]])
        >>> x.length_of_row(0)
        3
        >>> x.length_of_row(1)
        2
        """
        return len(self.row(row))

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
        >>> from landlab.utils.jaggedarray_ma import JaggedArray
        >>> x = JaggedArray([[0, 1, 2], [3, 4]])
        >>> x.row(0)
        array([0, 1, 2])
        >>> x.row(1)
        array([3, 4])
        """
        return self._values[row].compressed()

    def __iter__(self):
        """Iterate over the rows of the array.

        Examples
        --------
        >>> from landlab.utils.jaggedarray_ma import JaggedArray
        >>> x = JaggedArray([[0, 1, 2], [3, 4]])
        >>> for row in x: row
        array([0, 1, 2])
        array([3, 4])
        """
        for row in self._values:
            yield row.compressed()

    def foreach_row(self, func, out=None):
        """Apply an operator row-by-row.

        Examples
        --------
        >>> from landlab.utils.jaggedarray_ma import JaggedArray
        >>> x = JaggedArray([[0, 1, 2], [3, 4]])
        >>> x.foreach_row(np.sum)
        array([3, 7])

        >>> out = np.empty(2, dtype=int)
        >>> x.foreach_row(np.sum, out=out) is out
        True
        >>> out
        array([3, 7])
        """
        if out is None:
            return func(self._values, axis=1).compressed()
        else:
            return func(self._values, axis=1, out=out)

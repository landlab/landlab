#! /usr/bin/env python
"""Container that holds a collection of named data-fields."""

import numpy as np

_UNKNOWN_UNITS = "?"


class Error(Exception):

    """Base class for errors in this module."""

    pass


class FieldError(Error, KeyError):

    """Raise this error for a missing field name."""

    def __init__(self, field):
        self._field = field

    def __str__(self):
        return self._field


def need_to_reshape_array(array, field_size):
    """Check to see if an array needs to be resized before storing.

    When possible, a reference to an array is stored. However, if the
    array is not of the correct shape, a view of the array (with
    the correct shape) is stored.

    Parameters
    ----------
    array : numpy array
        Numpy array to check.
    field_size : int
        Size of the field the array will be placed into.

    Returns
    -------
    bool
        True is the array should be resized.
    """
    if field_size > 1:
        stored_shape = (field_size,)
    else:
        stored_shape = array.squeeze().shape

    return array.shape != stored_shape


class ScalarDataFields(dict):

    """Collection of named data fields that are of the same size.

    Holds a collection of data fields that all contain the same number of
    elements and index each of them with a name. This class inherits from
    a standard Python `dict`, which allows access to the fields through
    dict-like syntax.

    The syntax `.at_[element]` can also be used as syntactic sugar to access
    fields. e.g., `n1 = fields.at_node['name1']`, `n2 = grid.at_link['name2']`.

    Parameters
    ----------
    size : int
        The number of elements in each of the data fields.

    Attributes
    ----------
    units
    size

    See Also
    --------
    landlab.field.ModelDataFields.ones : Hold collections of
        `ScalarDataFields`.

    Examples
    --------
    >>> from landlab.field import ScalarDataFields
    >>> fields = ScalarDataFields(5)
    >>> fields.add_field('land_surface__elevation', [1, 2, 3, 4, 5])
    array([1, 2, 3, 4, 5])
    >>> fields['air__temperature'] = [2, 3, 4, 5, 6]
    >>> fields['land_surface__temperature'] = [0, 1]
    ...     # doctest: +IGNORE_EXCEPTION_DETAIL
    Traceback (most recent call last):
    ValueError: total size of the new array must be the same as the field

    Fields can also be multidimensional arrays so long as they can be
    resized such that the first dimension is the size of the field.
    The stored field will be resized view of the input array such that
    the size of the first dimension is the size of the field.

    >>> fields['air__temperature'] = [2, 3, 4, 5, 6, 7, 8, 9, 10, 11]
    >>> fields['air__temperature']
    array([[ 2,  3],
           [ 4,  5],
           [ 6,  7],
           [ 8,  9],
           [10, 11]])

    You can also create unsized fields. These fields will not be sized until
    the first field is added to the collection. Once the size is set, all
    fields must be the same size.

    >>> fields = ScalarDataFields()
    >>> fields['land_surface__temperature'] = [0, 1]
    >>> fields['land_surface__temperature']
    array([0, 1])
    >>> fields['air__temperature'] = [2, 3, 4, 5, 6]
    ...     # doctest: +IGNORE_EXCEPTION_DETAIL
    Traceback (most recent call last):
    ValueError: total size of the new array must be the same as the field

    Fields defined on a grid, which inherits from the ScalarModelFields class,
    behave similarly, though the length of those fields will be forced
    by the element type they are defined on:

    >>> from landlab import RasterModelGrid
    >>> import numpy as np
    >>> mg = RasterModelGrid((4, 5))
    >>> z = mg.add_field(
    ...     "topographic__elevation",
    ...     np.random.rand(mg.number_of_cells),
    ...     at="cell",
    ...     units="m",
    ... )
    >>> mg.at_cell['topographic__elevation'].size == mg.number_of_cells
    True

    LLCATS: FIELDCR, FIELDIO
    """

    def __init__(self, size=None):
        self._size = size

        super(ScalarDataFields, self).__init__()
        self._units = dict()

    @property
    def units(self):
        """Get units for values of the field.

        Returns
        -------
        str
            Units of the field.
        """
        return self._units

    @property
    def size(self):
        """Number of elements in the field.

        Returns
        -------
        int
            The number of elements in the field.
        """
        try:
            self._size
        except AttributeError:
            self._size = None
        return self._size

    @size.setter
    def size(self, size):
        if self._size is None:
            self._size = size
        else:
            raise ValueError("size has already been set")

    def empty(self, **kwds):
        """Uninitialized array whose size is that of the field.

        Return a new array of the data field size, without initializing
        entries. Keyword arguments are the same as that for the equivalent
        numpy function.

        See Also
        --------
        numpy.empty : See for a description of optional keywords.
        landlab.field.ScalarDataFields.ones : Equivalent method that
            initializes the data to 1.
        landlab.field.ScalarDataFields.zeros : Equivalent method that
            initializes the data to 0.

        Examples
        --------
        >>> from landlab.field import ScalarDataFields
        >>> field = ScalarDataFields(4)
        >>> field.empty() # doctest: +SKIP
        array([  2.31584178e+077,  -2.68156175e+154,   9.88131292e-324,
        ... 2.78134232e-309]) # Uninitialized memory

        Note that a new field is *not* added to the collection of fields.

        >>> list(field.keys())
        []
        """
        return np.empty(self.size, **kwds)

    def ones(self, **kwds):
        """Array, initialized to 1, whose size is that of the field.

        Return a new array of the data field size, filled with ones. Keyword
        arguments are the same as that for the equivalent numpy function.

        See Also
        --------
        numpy.ones : See for a description of optional keywords.
        landlab.field.ScalarDataFields.empty : Equivalent method that
            does not initialize the new array.
        landlab.field.ScalarDataFields.zeros : Equivalent method that
            initializes the data to 0.

        Examples
        --------
        >>> from landlab.field import ScalarDataFields
        >>> field = ScalarDataFields(4)
        >>> field.ones()
        array([ 1.,  1.,  1.,  1.])
        >>> field.ones(dtype=int)
        array([1, 1, 1, 1])

        Note that a new field is *not* added to the collection of fields.

        >>> list(field.keys())
        []
        """
        return np.ones(self.size, **kwds)

    def zeros(self, **kwds):
        """Array, initialized to 0, whose size is that of the field.

        Return a new array of the data field size, filled with zeros. Keyword
        arguments are the same as that for the equivalent numpy function.

        See Also
        --------
        numpy.zeros : See for a description of optional keywords.
        landlab.field.ScalarDataFields.empty : Equivalent method that does not
            initialize the new array.
        landlab.field.scalar_data_fields.ScalarDataFields.ones : Equivalent
            method that initializes the data to 1.

        Examples
        --------
        >>> from landlab.field import ScalarDataFields
        >>> field = ScalarDataFields(4)
        >>> field.zeros()
        array([ 0.,  0.,  0.,  0.])

        Note that a new field is *not* added to the collection of fields.

        >>> list(field.keys())
        []
        """
        return np.zeros(self.size, **kwds)

    def add_empty(self, name, units=_UNKNOWN_UNITS, clobber=False, **kwds):
        """Create and add an uninitialized array of values to the field.

        Create a new array of the data field size, without initializing
        entries, and add it to the field as *name*. The *units* keyword gives
        the units of the new fields as a string. Remaining keyword arguments
        are the same as that for the equivalent numpy function.

        Parameters
        ----------
        name : str
            Name of the new field to add.
        units : str, optional
            Optionally specify the units of the field.
        clobber : boolean, optional
            Raise an exception if adding to an already existing field.

        Returns
        -------
        array :
            A reference to the newly-created array.

        See Also
        --------
        numpy.empty : See for a description of optional keywords.
        landlab.field.ScalarDataFields.empty : Equivalent method that
            does not initialize the new array.
        landlab.field.ScalarDataFields.zeros : Equivalent method that
            initializes the data to 0.

        LLCATS: FIELDCR
        """
        return self.add_field(name, self.empty(**kwds), units=units, clobber=clobber)

    def add_ones(self, name, units=_UNKNOWN_UNITS, clobber=False, **kwds):
        """Create and add an array of values, initialized to 1, to the field.

        Create a new array of the data field size, filled with ones, and
        add it to the field as *name*. The *units* keyword gives the units of
        the new fields as a string. Remaining keyword arguments are the same
        as that for the equivalent numpy function.

        Parameters
        ----------
        name : str
            Name of the new field to add.
        units : str, optional
            Optionally specify the units of the field.
        clobber : boolean, optional
            Raise an exception if adding to an already existing field.

        Returns
        -------
        array :
            A reference to the newly-created array.

        See Also
        --------
        numpy.ones : See for a description of optional keywords.
        andlab.field.ScalarDataFields.add_empty : Equivalent method that
            does not initialize the new array.
        andlab.field.ScalarDataFields.add_zeros : Equivalent method that
            initializes the data to 0.

        Examples
        --------
        Add a new, named field to a collection of fields.

        >>> from landlab.field import ScalarDataFields
        >>> field = ScalarDataFields(4)
        >>> field.add_ones('topographic__elevation')
        array([ 1.,  1.,  1.,  1.])
        >>> list(field.keys())
        ['topographic__elevation']
        >>> field['topographic__elevation']
        array([ 1.,  1.,  1.,  1.])

        LLCATS: FIELDCR
        """
        return self.add_field(name, self.ones(**kwds), units=units, clobber=clobber)

    def add_zeros(self, name, units=_UNKNOWN_UNITS, clobber=False, **kwds):
        """Create and add an array of values, initialized to 0, to the field.

        Create a new array of the data field size, filled with zeros, and
        add it to the field as *name*. The *units* keyword gives the units of
        the new fields as a string. Remaining keyword arguments are the same
        as that for the equivalent numpy function.

        Parameters
        ----------
        name : str
            Name of the new field to add.
        units : str, optional
            Optionally specify the units of the field.
        clobber : boolean, optional
            Raise an exception if adding to an already existing field.

        Returns
        -------
        array :
            A reference to the newly-created array.

        See also
        --------
        numpy.zeros : See for a description of optional keywords.
        landlab.field.ScalarDataFields.add_empty : Equivalent method that
            does not initialize the new array.
        landlab.field.ScalarDataFields.add_ones : Equivalent method that
            initializes the data to 1.

        LLCATS: FIELDCR
        """
        return self.add_field(name, self.zeros(**kwds), units=units, clobber=clobber)

    def add_field(
        self, name, value_array, units=_UNKNOWN_UNITS, copy=False, clobber=False, **kwds
    ):
        """Add an array of values to the field.

        Add an array of data values to a collection of fields and associate it
        with the key, *name*. Use the *copy* keyword to, optionally, add a
        copy of the provided array.

        Parameters
        ----------
        name : str
            Name of the new field to add.
        value_array : numpy.array
            Array of values to add to the field.
        units : str, optional
            Optionally specify the units of the field.
        copy : boolean, optional
            If True, add a *copy* of the array to the field. Otherwise save add
            a reference to the array.
        clobber : boolean, optional
            Raise an exception if adding to an already existing field.

        Returns
        -------
        numpy.array
            The data array added to the field. Depending on the *copy*
            keyword, this could be a copy of *value_array* or *value_array*
            itself.

        Raises
        ------
        ValueError :
            If *value_array* has a size different from the field.

        Examples
        --------
        >>> import numpy as np
        >>> from landlab.field import ScalarDataFields
        >>> field = ScalarDataFields(4)
        >>> values = np.ones(4, dtype=int)
        >>> field.add_field('topographic__elevation', values)
        array([1, 1, 1, 1])

        A new field is added to the collection of fields. The saved value
        array is the same as the one initially created.

        >>> field['topographic__elevation'] is values
        True

        If you want to save a copy of the array, use the *copy* keyword. In
        addition, adding values to an existing field will remove the reference
        to the previously saved array. The *clobber=False* keyword changes this
        behavior to raise an exception in such a case.

        >>> field.add_field(
        ...     "topographic__elevation", values, copy=True, clobber=True
        ... )
        array([1, 1, 1, 1])
        >>> field['topographic__elevation'] is values
        False
        >>> field.add_field('topographic__elevation', values, clobber=False)
        ...     # doctest: +IGNORE_EXCEPTION_DETAIL
        Traceback (most recent call last):
        FieldError: topographic__elevation already exists

        LLCATS: FIELDCR
        """
        if not clobber and name in self:
            raise FieldError("{name}: already exists".format(name=name))

        value_array = np.asarray(value_array)

        if copy:
            value_array = value_array.copy()

        self[name] = value_array

        self.set_units(name, units)
        return self[name]

    def set_units(self, name, units):
        """Set the units for a field of values.

        Parameters
        ----------
        name: str
            Name of the field.
        units: str
            Units for the field

        Raises
        ------
        KeyError
            If the named field does not exist.

        LLCATS: FIELDCR
        """
        try:
            self._units[name] = units
        except AttributeError:
            self._units = {}
            self._units[name] = units

    def __setitem__(self, name, value_array):
        """Store a data field by name."""
        value_array = np.asarray(value_array)

        if self.size is None:
            self.size = value_array.size

        if need_to_reshape_array(value_array, self.size):
            value_array = value_array.reshape((self.size, -1)).squeeze()

        if name not in self:
            self.set_units(name, None)

        super(ScalarDataFields, self).__setitem__(name, value_array)

    def __getitem__(self, name):
        """Get a data field by name."""
        try:
            return super(ScalarDataFields, self).__getitem__(name)
        except KeyError:
            raise FieldError(name)

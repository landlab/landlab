#! /usr/bin/env python

import numpy as np


_UNKNOWN_UNITS = '?'


class ScalarDataFields(dict):
    def __init__(self, size):
        self._size = size

        super(ScalarDataFields, self).__init__()
        self._units = dict()

    @property
    def units(self):
        """Units for values of the field.

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
        return self._size

    def empty(self, **kwds):
        """
        Return a new array of the data field size, without initializing
        entries. Keyword arguments are the same as that for the equivalent
        numpy function.
        """
        return np.empty(self.size, **kwds)

    def ones(self, **kwds):
        """
        Return a new array of the data field size, filled with ones. Keyword
        arguments are the same as that for the equivalent numpy function.
        """
        return np.ones(self.size, **kwds)

    def zeros(self, **kwds):
        """
        Return a new array of the data field size, filled with zeros. Keyword
        arguments are the same as that for the equivalent numpy function.
        """
        return np.zeros(self.size, **kwds)

    def add_empty(self, name, units=_UNKNOWN_UNITS, **kwds):
        """Create and add an uninitialized array of values to the field.

        Parameters
        ----------
        name : str
            Name of the new field to add.
        units : str, optional
            Optionally specify the units of the field.

        Create a new array of the data field size, without initializing
        entries, and add it to the field as *name*. The *units* keyword gives
        the units of the new fields as a string. Remaining keyword arguments
        are the same as that for the equivalent numpy function.

        Returns a reference to the newly-created array.
        """
        return self.add_field(name, self.empty(**kwds), units=units)

    def add_ones(self, name, units=_UNKNOWN_UNITS, **kwds):
        """Create and add an array of values, initialized to 1, to the field.

        Parameters
        ----------
        name : str
            Name of the new field to add.
        units : str, optional
            Optionally specify the units of the field.

        Create a new array of the data field size, filled with ones, and
        add it to the field as *name*. The *units* keyword gives the units of
        the new fields as a string. Remaining keyword arguments are the same
        as that for the equivalent numpy function.

        Returns a reference to the newly-created array.
        """
        return self.add_field(name, self.ones(**kwds), units=units)

    def add_zeros(self, name, units=_UNKNOWN_UNITS, **kwds):
        """Create and add an array of values, initialized to 0, to the field.

        Parameters
        ----------
        name : str
            Name of the new field to add.
        units : str, optional
            Optionally specify the units of the field.

        Create a new array of the data field size, filled with zeros, and
        add it to the field as *name*. The *units* keyword gives the units of
        the new fields as a string. Remaining keyword arguments are the same
        as that for the equivalent numpy function.

        Returns a reference to the newly-created array.
        """
        return self.add_field(name, self.zeros(**kwds), units=units)

    def add_field(self, name, value_array, units=_UNKNOWN_UNITS, copy=False):
        """Add an array of values to the field.

        Parameters
        ----------
        name : str
            Name of the new field to add.
        value_array : numpy.array
            Array of values to add to the field.
        units : str, optional
            Optionally specify the units of the field.

        Returns
        -------
        numpy.array
            The input *value_array*.
        """
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
        """
        self._units[name] = units

    def __setitem__(self, name, value_array):
        assert(value_array.size == self.size)

        if name not in self:
            self.set_units(name, None)

        super(ScalarDataFields, self).__setitem__(name, value_array)

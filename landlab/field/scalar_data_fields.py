#! /usr/bin/env python

import numpy as np


_UNKNOWN_UNITS = '?'


class ScalarDataFields(dict):
    def __init__(self, size):
        self._size = size

        #print 'ScalarDataFields.__init__'
        super(ScalarDataFields, self).__init__()
        self._units = dict()

    @property
    def units(self):
        return self._units

    @property
    def size(self):
        return self._size

    @property
    def number_of_values(self):
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
        """
        Create a new array of the data field size, without initializing
        entries, and add it to the field as *name*. The *units* keyword gives
        the units of the new fields as a string. Remaining keyword arguments
        are the same as that for the equivalent numpy function.

        Returns a reference to the newly-created array.
        """
        return self.add_field(name, self.empty(**kwds), units=units)

    def add_ones(self, name, units=_UNKNOWN_UNITS, **kwds):
        """
        Create a new array of the data field size, filled with ones, and
        add it to the field as *name*. The *units* keyword gives the units of
        the new fields as a string. Remaining keyword arguments are the same
        as that for the equivalent numpy function.

        Returns a reference to the newly-created array.
        """
        return self.add_field(name, self.ones(**kwds), units=units)

    def add_zeros(self, name, units=_UNKNOWN_UNITS, **kwds):
        """
        Create a new array of the data field size, filled with zeros, and
        add it to the field as *name*. The *units* keyword gives the units of
        the new fields as a string. Remaining keyword arguments are the same
        as that for the equivalent numpy function.

        Returns a reference to the newly-created array.
        """
        return self.add_field(name, self.zeros(**kwds), units=units)

    def add_field(self, name, value_array, units=_UNKNOWN_UNITS, copy=False):
        if copy:
            value_array = value_array.copy()

        self[name] = value_array

        self.set_units(name, units)
        return self[name]

    def set_units(self, name, units):
        self._units[name] = units

    def __setitem__(self, name, value_array):
        assert(value_array.size == self.size)

        if name not in self:
            self.set_units(name, None)

        super(ScalarDataFields, self).__setitem__(name, value_array)

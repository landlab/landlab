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
        return self._units

    @property
    def size(self):
        return self._size

    def empty(self, **kwds):
        return np.empty(self.size, **kwds)

    def ones(self, **kwds):
        return np.ones(self.size, **kwds)

    def zeros(self, **kwds):
        return np.zeros(self.size, **kwds)

    def add_empty(self, name, units=_UNKNOWN_UNITS, **kwds):
        self.add_field(name, self.empty(**kwds), units=units)

    def add_ones(self, name, units=_UNKNOWN_UNITS, **kwds):
        self.add_field(name, self.ones(**kwds), units=units)

    def add_zeros(self, name, units=_UNKNOWN_UNITS, **kwds):
        self.add_field(name, self.zeros(**kwds), units=units)

    def add_field(self, name, value_array, units=_UNKNOWN_UNITS, copy=False):
        if copy:
            value_array = value_array.copy()

        self[name] = value_array

        self.set_units(name, units)

    def set_units(self, name, units):
        self._units[name] = units

    def __setitem__(self, name, value_array):
        assert(value_array.size == self.size)

        if name not in self:
            self.set_units(name, None)

        super(ScalarDataFields, self).__setitem__(name, value_array)

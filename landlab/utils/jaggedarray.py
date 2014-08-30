import numpy as np


class JaggedArray(object):
    def __init__(self, values, values_per_row):
        self._values = values
        self._number_of_rows = len(values_per_row)
        self._offsets = JaggedArray._offsets_from_values_per_row(values_per_row)

    @property
    def size(self):
        return len(self._values)

    @property
    def number_of_rows(self):
        return self._number_of_rows

    @staticmethod
    def _offsets_from_values_per_row(values_per_row):
        offset = np.empty(len(values_per_row) + 1, dtype=int)
        np.cumsum(values_per_row, out=offset[1:])
        offset[0] = 0
        return offset

    def length_of_row(self, row):
        return self._offsets[row + 1] - self._offsets[row]

    def row(self, row):
        return self._values[self._offsets[row]:self._offsets[row + 1]]

    def iter(self):
        for n in xrange(self._number_of_rows):
            yield self.row(n)


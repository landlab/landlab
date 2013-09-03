#! /usr/bin/env python

import numpy as np
from landlab.field import ModelFields

class ModelStructuredFields(ModelFields):
    """
    >>> import numpy as np
    >>> field = ModelStructuredFields((4, 3))
    >>> values = np.arange(12.)
    >>> field.add_field('surface__elevation', values)
    >>> field.keys()
    ['surface__elevation']
    >>> print field['surface__elevation']
    [  0.   1.   2.   3.   4.   5.   6.   7.   8.   9.  10.  11.]
    >>> field.add_empty('new_var', dtype=np.int)
    """
    def __init__(self, shape):
        assert(len(shape) == 2)

        self._shape = tuple(shape)
        super(ModelStructuredFields, self).__init__(np.prod(shape))

    @property
    def shape(self):
        return self._shape

    def add_field(self, name, value_array, reshape=False, units=None,
                  copy=False):
        if copy:
            value_array = value_array.copy()

        if reshape:
            value_array = value_array.view()
            value_array.shape = self.shape

        super(ModelStructuredFields, self).add_field(name, value_array)

    def imshow(self, name, **kwds):
        from landlab.plot.imshow import imshow_field
        imshow_field(self, name, **kwds)

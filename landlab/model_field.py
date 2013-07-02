#! /usr/bin/env python

from landlab import RasterModelGrid


class RasterModelField(RasterModelGrid, dict):
    """
    >>> import numpy as np
    >>> field = RasterModelField((4, 3))
    >>> values = np.arange(12.)
    >>> field.add_field('surface__elevation', values)
    >>> field.keys()
    ['surface__elevation']
    """
    def __init__(self, shape, *args):
        try:
            origin = args[1]
        except IndexError:
            origin = (0., 0.)

        try:
            spacing = args[0]
        except IndexError:
            spacing = (1., 1.)

        assert(len(shape) == 2)
        assert(len(spacing) == 2)
        assert(len(origin) == 2)

        super(RasterModelField, self).__init__(num_rows=shape[0],
                                               num_cols=shape[1],
                                               dx=spacing[0])

    def add_field(self, name, values, centering='point', reshape=False):
        assert(values.size == self.get_count_of_all_cells())

        self[name] = values
        if reshape:
            self[name].shape = (self.get_count_of_rows(),
                                self.get_count_of_cols())

    def calculate_gradients_at_active_links(self, var_name):
        grad_name = var_name + '_gradient'

        try:
            super(RasterModelField, self).calculate_gradients_at_active_links(
                self[var_name], gradient=self[grad_name])
        except KeyError:
            self.add_field[grad_name] = (
                super(RasterModelField, self).calculate_gradients_at_active_links(
                    self[var_name]))

        return self[grad_name]


if __name__ == '__main__':
    import doctest
    doctest.testmod()

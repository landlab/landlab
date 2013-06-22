#!/usr/bin/env python

import numpy as np

from landlab.model_field import RasterModelField
from landlab.components.flexure.funcs import subside_grid


class Flexure(RasterModelField):
    """
    Landlab component that implements a 1 and 2D lithospheric flexure
    model.

    >>> flex = Flexure((5, 4), (1.e4, 1.e4), (0., 0.))
    >>> flex.name
    'Flexure'
    >>> sorted(flex.input_var_names)
    ['lithosphere__elevation', 'lithosphere__overlying_pressure']
    >>> sorted(flex.output_var_names)
    ['lithosphere__elevation', 'lithosphere__elevation_increment']
    >>> for var in sorted(flex.units): flex.units[var]
    'm'
    'm'
    'Pa'

    >>> flex.get_count_of_rows()
    5
    >>> flex.get_count_of_cols()
    4

    >>> np.all(flex['lithosphere__elevation'] == 0.)
    True
    >>> np.all(flex['lithosphere__overlying_pressure'] == 0.)
    True
    >>> flex.update()
    >>> np.all(flex['lithosphere__elevation_increment'] == 0.)
    True

    >>> load = flex['lithosphere__overlying_pressure']
    >>> load[2, 2] = 1e9
    >>> dz = flex['lithosphere__elevation_increment']
    >>> np.all(dz == 0.)
    True

    >>> flex.update()
    >>> np.all(flex['lithosphere__elevation_increment'] == 0.)
    False
    """

    _name = 'Flexure'

    _input_var_names = [
        'lithosphere__overlying_pressure',
        'lithosphere__elevation',
    ]
    _output_var_names = [
        'lithosphere__elevation_increment',
        'lithosphere__elevation',
    ]

    _var_units = {
        'lithosphere__overlying_pressure': 'Pa',
        'lithosphere__elevation': 'm',
        'lithosphere__elevation_increment': 'm',
    }

    def __init__(self, shape, spacing, origin):
        super(Flexure, self).__init__(shape, spacing)

        self._shape = shape
        self._spacing = spacing
        self._origin = origin

        for name in self._input_var_names:
            self.add_field(name, np.zeros(shape[0] * shape[1], dtype=np.float))

        for name in self._output_var_names:
            self.add_field(name, np.zeros(shape[0] * shape[1], dtype=np.float))

        self._last_load = self['lithosphere__overlying_pressure'].copy()
        (self._x, self._y) = np.meshgrid(
            np.linspace(origin[0], origin[0] + (shape[0] + 1) * spacing[0]),
            np.linspace(origin[1], origin[1] + (shape[1] + 1) * spacing[1]))

        self._eet = 65000.
        self._youngs = 7e10

    @property
    def input_var_names(self):
        return self._input_var_names

    @property
    def output_var_names(self):
        return self._output_var_names

    @property
    def name(self):
        return self._name

    @property
    def units(self):
        return self._var_units

    #def initialize(self):
    #    pass

    def update(self):
        elevation = self['lithosphere__elevation']
        load = self['lithosphere__overlying_pressure']
        deflection = self['lithosphere__elevation_increment']

        new_load = load - self._last_load

        self._last_load = load.copy()

        (x_coords, y_coords) = (self.get_cell_x_coords(),
                                self.get_cell_y_coords())

        x_coords.shape = deflection.shape
        y_coords.shape = deflection.shape

        subside_grid(deflection, x_coords, y_coords, new_load, self._eet,
                     self._youngs)

        elevation += deflection

    def finalize(self):
        pass


if __name__ == "__main__":
    import doctest
    doctest.testmod()

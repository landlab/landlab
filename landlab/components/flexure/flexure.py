#!/usr/bin/env python

import numpy as np

from landlab.model_field import RasterModelField
from landlab.components.flexure.funcs import subside_point_loads


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
    >>> load[4] = 1e9
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

    @property
    def shape(self):
        return self._shape

    @property
    def coords(self):
        return (self.get_cell_x_coords(), self.get_cell_y_coords())

    def update(self, n_procs=1):
        elevation = self['lithosphere__elevation']
        load = self['lithosphere__overlying_pressure']
        deflection = self['lithosphere__elevation_increment']

        new_load = load - self._last_load

        self._last_load = load.copy()

        deflection.fill(0.)
        self.subside_loads(new_load, self.coords, deflection=deflection,
                           n_procs=n_procs)

        elevation += deflection

    def subside_loads(self, loads, locs, deflection=None, n_procs=1):
        if deflection is None:
            deflection = np.empty(self.shape, dtype=np.float)

        subside_point_loads(loads, locs, self.coords, self._eet, self._youngs,
                            deflection=deflection, n_procs=n_procs)

        return deflection

    def subside_load(self, load, loc, deflection=None):
        subside_point_load(load, loc, self.coords, self._eet, self._youngs,
                           deflection=self['lithosphere__elevation_increment'])

        return deflection


if __name__ == "__main__":
    import doctest
    doctest.testmod()

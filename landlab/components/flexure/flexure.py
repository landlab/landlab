#!/usr/bin/env python
"""Deform the lithosphere with 1D or 2D flexure.

Landlab component that implements a 1 and 2D lithospheric flexure
model.

Examples
--------

Create a grid on which we will run the flexure calculations.

>>> from landlab import RasterModelGrid
>>> from landlab.components.flexure import Flexure
>>> grid = RasterModelGrid((5, 4), spacing=(1.e4, 1.e4))

Check the fields that are used as input to the flexure component.

>>> Flexure.input_var_names # doctest: +NORMALIZE_WHITESPACE
('lithosphere__overlying_pressure', 'lithosphere__elevation',
 'planet_surface_sediment__deposition_increment')

Check the units for the fields.

>>> Flexure.var_units('lithosphere__elevation')
'm'
>>> Flexure.var_units('lithosphere__overlying_pressure')
'Pa'

If you are not sure about one of the input or output variables, you can
get help for specific variables.

>>> Flexure.var_help('planet_surface_sediment__deposition_increment')
name: planet_surface_sediment__deposition_increment
description:
  The amount of sediment deposited at the land surface in one timestep
units: m
at: node
intent: in

Add these fields to the grid.

>>> _ = grid.add_zeros('lithosphere__elevation', at='node')
>>> _ = grid.add_zeros('lithosphere__overlying_pressure', at='node')
>>> _ = grid.add_zeros('planet_surface_sediment__deposition_increment', at='node')

>>> dh = grid.at_node['planet_surface_sediment__deposition_increment']
>>> dh = dh.reshape(grid.shape)
>>> dh[1:-1, 1:-1] = 3300 / 2650.

>>> flex = Flexure(grid)
>>> flex.update()

>>> flex.output_var_names
('lithosphere__elevation_increment', 'lithosphere__elevation')
>>> flex.grid.at_node['lithosphere__elevation']
...     # doctest: +NORMALIZE_WHITESPACE
array([ 0.,  0.,  0.,  0.,
        0., -1., -1.,  0.,
        0., -1., -1.,  0.,
        0., -1., -1.,  0.,
        0.,  0.,  0.,  0.])
"""

import numpy as np

from landlab import Component
from .funcs import get_flexure_parameter
from ...utils.decorators import use_file_name_or_kwds


class Flexure(Component):

    """Deform the lithosphere with 1D or 2D flexure.

    Landlab component that implements a 1 and 2D lithospheric flexure
    model.

    Examples
    --------
    >>> from landlab import RasterModelGrid
    >>> from landlab.components.flexure import Flexure
    >>> grid = RasterModelGrid((5, 4), spacing=(1.e4, 1.e4))

    >>> flex = Flexure(grid)
    >>> flex.name
    'Flexure'
    >>> sorted(flex.input_var_names) # doctest: +NORMALIZE_WHITESPACE
    ['lithosphere__elevation', 'lithosphere__overlying_pressure',
     'planet_surface_sediment__deposition_increment']
    >>> sorted(flex.output_var_names)
    ['lithosphere__elevation', 'lithosphere__elevation_increment']
    >>> sorted(flex.units) # doctest: +NORMALIZE_WHITESPACE
    [('lithosphere__elevation', 'm'),
     ('lithosphere__elevation_increment', 'm'),
     ('lithosphere__overlying_pressure', 'Pa'),
     ('planet_surface_sediment__deposition_increment', 'm')]

    >>> flex.grid.number_of_node_rows
    5
    >>> flex.grid.number_of_node_columns
    4
    >>> flex.grid is grid
    True

    >>> np.all(grid.at_node['lithosphere__elevation'] == 0.)
    True

    >>> np.all(grid.at_node['lithosphere__elevation'] == 0.)
    True
    >>> np.all(grid.at_node['lithosphere__overlying_pressure'] == 0.)
    True
    >>> flex.update()
    >>> np.all(grid.at_node['lithosphere__elevation_increment'] == 0.)
    True

    >>> load = grid.at_node['lithosphere__overlying_pressure']
    >>> load[4] = 1e9
    >>> dz = grid.at_node['lithosphere__elevation_increment']
    >>> np.all(dz == 0.)
    True

    >>> flex.update()
    >>> np.all(grid.at_node['lithosphere__elevation_increment'] == 0.)
    False
    """

    _name = 'Flexure'

    _input_var_names = (
        'lithosphere__overlying_pressure',
        'lithosphere__elevation',
        'planet_surface_sediment__deposition_increment',
    )

    _output_var_names = (
        'lithosphere__elevation_increment',
        'lithosphere__elevation',
    )

    _var_units = {
        'lithosphere__overlying_pressure': 'Pa',
        'lithosphere__elevation': 'm',
        'lithosphere__elevation_increment': 'm',
        'planet_surface_sediment__deposition_increment': 'm',
    }

    _var_mapping = {
        'lithosphere__overlying_pressure': 'node',
        'lithosphere__elevation': 'node',
        'lithosphere__elevation_increment': 'node',
        'planet_surface_sediment__deposition_increment': 'node',
    }

    _var_doc = {
        'lithosphere__overlying_pressure':
            'The pressure at the base of the lithosphere',
        'lithosphere__elevation':
            'The elevation of the top of the lithosphere, i.e., the land '
            'surface',
        'lithosphere__elevation_increment':
            'The change in elevation of the top of the lithosphere (the land '
            'surface) in one timestep',
        'planet_surface_sediment__deposition_increment':
            'The amount of sediment deposited at the land surface in one '
            'timestep',
    }

    @use_file_name_or_kwds
    def __init__(self, grid, eet=65e3, youngs=7e10, method='airy',
                 rho_sed=2650., rho_mantle=3300., gravity=9.80665, **kwds):
        """Initialize the flexure component.

        Parameters
        ----------
        grid : RasterModelGrid
            A grid.
        eet : float, optional
            Effective elastic thickness (m).
        youngs : float, optional
            Young's modulus.
        method : {'airy', 'flexure'}, optional
            Method to use to calculate deflections.
        """
        if method not in ('airy', 'flexure'):
            raise ValueError(
                '{method}: method not understood'.format(method=method))

        self._eet = eet
        self._youngs = youngs
        self._method = method
        self._rho_sed = rho_sed
        self._rho_mantle = rho_mantle
        self._gravity = gravity

        self._grid = grid

        super(Flexure, self).__init__(grid, **kwds)

        for name in self._input_var_names:
            if name not in self.grid.at_node:
                self.grid.add_zeros('node', name, units=self._var_units[name])

        for name in self._output_var_names:
            if name not in self.grid.at_node:
                self.grid.add_zeros('node', name, units=self._var_units[name])

        self._last_load = (
            self.grid.at_node['lithosphere__overlying_pressure'].copy())

        self._r = self._set_kei_func_grid()

    @property
    def eet(self):
        return self._eet

    @property
    def youngs(self):
        return self._youngs

    @property
    def rho_sed(self):
        return self._rho_sed

    @property
    def gamma_sed(self):
        return self._rho_sed * self._gravity

    @property
    def rho_mantle(self):
        return self._rho_mantle

    @property
    def gamma_mantle(self):
        return self._rho_mantle * self._gravity

    @property
    def gravity(self):
        return self._gravity

    @property
    def method(self):
        return self._method

    @property
    def alpha(self):
        return get_flexure_parameter(self._eet, self._youngs, 2,
                                     gamma_mantle=self.gamma_mantle)

    def _set_kei_func_grid(self):
        from scipy.special import kei

        dx, dy = np.meshgrid(
            np.arange(self._grid.number_of_node_columns) * self._grid.dx,
            np.arange(self._grid.number_of_node_rows) * self._grid.dy)

        return kei(np.sqrt(dx ** 2 + dy ** 2) / self.alpha)

    def update(self, n_procs=1):
        elevation = self.grid.at_node['lithosphere__elevation']
        load = self.grid.at_node['lithosphere__overlying_pressure']
        deflection = self.grid.at_node['lithosphere__elevation_increment']
        deposition = self.grid.at_node['planet_surface_sediment__deposition_increment']

        new_load = ((load - self._last_load) +
                    (deposition * self.gamma_sed).flat)

        self._last_load = load.copy()

        deflection.fill(0.)

        if self._method == 'airy':
            deflection[:] = new_load / self.gamma_mantle
        else:
            self.subside_loads(new_load, deflection=deflection,
                               n_procs=n_procs)
        np.subtract(elevation, deflection, out=elevation)

    def subside_loads(self, loads, deflection=None, n_procs=1):
        if deflection is None:
            deflection = np.empty(self.shape, dtype=np.float)

        from .cfuncs import subside_grid_in_parallel

        w = deflection.reshape(self._grid.shape)
        load = loads.reshape(self._grid.shape)

        subside_grid_in_parallel(w, load, self._r, self.alpha,
                                 self.gamma_mantle, n_procs)

        return deflection

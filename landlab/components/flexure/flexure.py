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
('lithosphere__overlying_pressure_increment',)

Check the units for the fields.

>>> Flexure.var_units('lithosphere__overlying_pressure_increment')
'Pa'

If you are not sure about one of the input or output variables, you can
get help for specific variables.

>>> Flexure.var_help('lithosphere__overlying_pressure_increment')
name: lithosphere__overlying_pressure_increment
description:
  Applied pressure to the lithosphere over a time step
units: Pa
at: node
intent: in

>>> flex = Flexure(grid)

In creating the component, a field (initialized with zeros) was added to the
grid. Reset the interior nodes for the loading.

>>> dh = grid.at_node['lithosphere__overlying_pressure_increment']
>>> dh = dh.reshape(grid.shape)
>>> dh[1:-1, 1:-1] = flex.gamma_mantle

>>> flex.update()

>>> flex.output_var_names
('lithosphere_surface__elevation_increment',)
>>> flex.grid.at_node['lithosphere_surface__elevation_increment']
...     # doctest: +NORMALIZE_WHITESPACE
array([ 0., 0., 0., 0.,
        0., 1., 1., 0.,
        0., 1., 1., 0.,
        0., 1., 1., 0.,
        0., 0., 0., 0.])
"""

import numpy as np

from landlab import Component
from .funcs import get_flexure_parameter
from ...utils.decorators import use_file_name_or_kwds


class Flexure(Component):

    """Deform the lithosphere with 1D or 2D flexure.

    Landlab component that implements a 1 and 2D lithospheric flexure
    model.

    Construction::

        Flexure(grid, eet=65e3, youngs=7e10, method='airy', rho_mantle=3300.,
                gravity=9.80665)

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
    rho_mantle : float, optional
        Density of the mantle (kg / m^3).
    gravity : float, optional
        Acceleration due to gravity (m / s^2).

    Examples
    --------
    >>> from landlab import RasterModelGrid
    >>> from landlab.components.flexure import Flexure
    >>> grid = RasterModelGrid((5, 4), spacing=(1.e4, 1.e4))

    >>> flex = Flexure(grid)
    >>> flex.name
    'Flexure'
    >>> flex.input_var_names
    ('lithosphere__overlying_pressure_increment',)
    >>> flex.output_var_names
    ('lithosphere_surface__elevation_increment',)
    >>> sorted(flex.units) # doctest: +NORMALIZE_WHITESPACE
    [('lithosphere__overlying_pressure_increment', 'Pa'),
     ('lithosphere_surface__elevation_increment', 'm')]

    >>> flex.grid.number_of_node_rows
    5
    >>> flex.grid.number_of_node_columns
    4
    >>> flex.grid is grid
    True

    >>> np.all(grid.at_node['lithosphere_surface__elevation_increment'] == 0.)
    True

    >>> np.all(grid.at_node['lithosphere__overlying_pressure_increment'] == 0.)
    True
    >>> flex.update()
    >>> np.all(grid.at_node['lithosphere_surface__elevation_increment'] == 0.)
    True

    >>> load = grid.at_node['lithosphere__overlying_pressure_increment']
    >>> load[4] = 1e9
    >>> dz = grid.at_node['lithosphere_surface__elevation_increment']
    >>> np.all(dz == 0.)
    True

    >>> flex.update()
    >>> np.all(grid.at_node['lithosphere_surface__elevation_increment'] == 0.)
    False
    """

    _name = 'Flexure'

    _input_var_names = (
        'lithosphere__overlying_pressure_increment',
    )

    _output_var_names = (
        'lithosphere_surface__elevation_increment',
    )

    _var_units = {
        'lithosphere__overlying_pressure_increment': 'Pa',
        'lithosphere_surface__elevation_increment': 'm',
    }

    _var_mapping = {
        'lithosphere__overlying_pressure_increment': 'node',
        'lithosphere_surface__elevation_increment': 'node',
    }

    _var_doc = {
        'lithosphere__overlying_pressure_increment':
            'Applied pressure to the lithosphere over a time step',
        'lithosphere_surface__elevation_increment':
            'The change in elevation of the top of the lithosphere (the land '
            'surface) in one timestep',
    }

    @use_file_name_or_kwds
    def __init__(self, grid, eet=65e3, youngs=7e10, method='airy',
                 rho_mantle=3300., gravity=9.80665, **kwds):
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
        rho_mantle : float, optional
            Density of the mantle (kg / m^3).
        gravity : float, optional
            Acceleration due to gravity (m / s^2).
        """
        if method not in ('airy', 'flexure'):
            raise ValueError(
                '{method}: method not understood'.format(method=method))

        self._eet = eet
        self._youngs = youngs
        self._method = method
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

        self._r = self._set_kei_func_grid()

    @property
    def eet(self):
        """Effective elastic thickness (m)."""
        return self._eet

    @eet.setter
    def eet(self, new_val):
        if new_val <= 0:
            raise ValueError('Effective elastic thickness must be positive.')
        self._r = self._set_kei_func_grid()
        self._eet = new_val

    @property
    def youngs(self):
        """Young's modulus of lithosphere (Pa)."""
        return self._youngs

    @property
    def rho_mantle(self):
        """Density of mantle (kg/m^3)."""
        return self._rho_mantle

    @property
    def gamma_mantle(self):
        """Specific density of mantle (N/m^3)."""
        return self._rho_mantle * self._gravity

    @property
    def gravity(self):
        """Acceleration due to gravity (m/s^2)."""
        return self._gravity

    @property
    def method(self):
        """Name of method used to calculate deflections."""
        return self._method

    @property
    def alpha(self):
        """Flexure parameter (m)."""
        return get_flexure_parameter(self._eet, self._youngs, 2,
                                     gamma_mantle=self.gamma_mantle)

    def _set_kei_func_grid(self):
        from scipy.special import kei

        dx, dy = np.meshgrid(
            np.arange(self._grid.number_of_node_columns) * self._grid.dx,
            np.arange(self._grid.number_of_node_rows) * self._grid.dy)

        return kei(np.sqrt(dx ** 2 + dy ** 2) / self.alpha)

    def update(self, n_procs=1):
        """Update fields with current loading conditions.

        Parameters
        ----------
        n_procs : int, optional
            Number of processors to use for calculations.
        """
        load = self.grid.at_node['lithosphere__overlying_pressure_increment']
        deflection = self.grid.at_node['lithosphere_surface__elevation_increment']

        new_load = load.copy()

        deflection.fill(0.)

        if self._method == 'airy':
            deflection[:] = new_load / self.gamma_mantle
        else:
            self.subside_loads(new_load, deflection=deflection,
                               n_procs=n_procs)

    def subside_loads(self, loads, deflection=None, n_procs=1):
        """Subside surface due to multiple loads.

        Parameters
        ----------
        loads : ndarray of float
            Loads applied to each grid node.
        deflection : ndarray of float, optional
            Buffer to place resulting deflection values.
        n_procs : int, optional
            Number of processors to use for calculations.

        Returns
        -------
        ndarray of float
            Deflections caused by the loading.
        """
        if deflection is None:
            deflection = np.empty(self.shape, dtype=np.float)

        from .cfuncs import subside_grid_in_parallel

        w = deflection.reshape(self._grid.shape)
        load = loads.reshape(self._grid.shape)

        subside_grid_in_parallel(w, load * self._grid.dx * self._grid.dy,
                                 self._r, self.alpha, self.gamma_mantle,
                                 n_procs)

        return deflection

#!/usr/bin/env python
"""Deform the lithosphere with 1D or 2D flexure.

Landlab component that implements a 1D lithospheric flexure.

Examples
--------

Create a grid on which we will run the flexure calculations.

>>> from landlab import RasterModelGrid
>>> from landlab.components.flexure import Flexure1D
>>> grid = RasterModelGrid((3, 4), xy_spacing=(1.0e4, 1.0e4))
>>> lith_press = grid.add_zeros("node", "lithosphere__increment_of_overlying_pressure")

Because `Flexure1D` is a one-dimensional component, it operates
*independently* on each row of grid nodes. By default, it will
calculate deflections on every row but you can also specify
which rows to operate on by passing an indexing array to the
`rows` keyword.

In this example, we'll just calculate deflections along the
middle row or nodes (that is, row 1)

>>> flex = Flexure1D(grid, rows=1)

In creating the component, a field (initialized with zeros) was
added to the grid that represents the current loading distribution.
If the grid already had this field, the component would use the
existing field. This can be accessed either through the *grid*
attribute in the same way as other landlab fields,

>>> flex.grid.at_node["lithosphere__increment_of_overlying_pressure"]
array([0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.])

or through the `load_at_node` attribute of `Flexure1D`,

>>> flex.load_at_node
array([[0., 0., 0., 0.],
       [0., 0., 0., 0.],
       [0., 0., 0., 0.]])

Notice that `load_at_node` returns a reshaped view of the array
whereas the field returns a flattened array. Change values in this
array to add loads to the grid,

>>> flex.load_at_node[1, 2:4] = flex.gamma_mantle
>>> flex.run_one_step()

The output deflections can be retrieved either using landlab fields
as,

>>> flex.grid.at_node["lithosphere_surface__increment_of_elevation"]
array([0., 0., 0., 0., 0., 0., 1., 1., 0., 0., 0., 0.])

or through the `dz_at_node` attribute,

>>> flex.dz_at_node
array([[0., 0., 0., 0.],
       [0., 0., 1., 1.],
       [0., 0., 0., 0.]])
"""
import contextlib

import numpy as np

from landlab import Component
from landlab.components.flexure._ext.flexure1d import subside_load_1d


class Flexure1D(Component):
    """Deform the lithosphere with 1D flexure.

    Landlab component that implements a 1D lithospheric flexure model.

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
    rho_water : float, optional
        Density of the sea water (kg / m^3).
    gravity : float, optional
        Acceleration due to gravity (m / s^2).
    rows : int, optional
        Node rows that this component will operate on (default is to
        operate on *all* rows).

    Examples
    --------
    >>> from landlab import RasterModelGrid
    >>> from landlab.components.flexure import Flexure1D
    >>> grid = RasterModelGrid((5, 4), xy_spacing=(1.0e4, 1.0e4))
    >>> lith_press = grid.add_zeros(
    ...     "node", "lithosphere__increment_of_overlying_pressure"
    ... )
    >>> flex = Flexure1D(grid)
    >>> flex.name
    '1D Flexure Equation'
    >>> flex.input_var_names
    ('lithosphere__increment_of_overlying_pressure',)
    >>> flex.output_var_names
    ('lithosphere_surface__increment_of_elevation',)
    >>> sorted(flex.units)
    [('lithosphere__increment_of_overlying_pressure', 'Pa'),
     ('lithosphere_surface__increment_of_elevation', 'm')]

    >>> flex.grid.number_of_node_rows
    5
    >>> flex.grid.number_of_node_columns
    4
    >>> flex.grid is grid
    True

    >>> np.all(grid.at_node["lithosphere_surface__increment_of_elevation"] == 0.0)
    True

    >>> np.all(grid.at_node["lithosphere__increment_of_overlying_pressure"] == 0.0)
    True
    >>> flex.update()
    >>> np.all(grid.at_node["lithosphere_surface__increment_of_elevation"] == 0.0)
    True

    >>> load = grid.at_node["lithosphere__increment_of_overlying_pressure"]
    >>> load[4] = 1e9
    >>> dz = grid.at_node["lithosphere_surface__increment_of_elevation"]
    >>> np.all(dz == 0.0)
    True

    >>> flex.update()
    >>> np.all(grid.at_node["lithosphere_surface__increment_of_elevation"] == 0.0)
    False
    """

    _name = "1D Flexure Equation"

    _unit_agnostic = True

    _info = {
        "lithosphere__increment_of_overlying_pressure": {
            "dtype": float,
            "intent": "in",
            "optional": False,
            "units": "Pa",
            "mapping": "node",
            "doc": "Applied pressure to the lithosphere over a time step",
        },
        "lithosphere_surface__increment_of_elevation": {
            "dtype": float,
            "intent": "out",
            "optional": False,
            "units": "m",
            "mapping": "node",
            "doc": (
                "The change in elevation of the top of the lithosphere (the "
                "land surface) in one timestep"
            ),
        },
    }

    _POISSON = 0.25

    def __init__(
        self,
        grid,
        eet=65e3,
        youngs=7e10,
        method="airy",
        rho_mantle=3300.0,
        rho_water=1030.0,
        gravity=9.80665,
        rows=None,
    ):
        """Initialize the flexure component.

        Parameters
        ----------
        grid : RasterModelGrid
            A grid.
        eet : float, optional
            Effective elastic thickness (m).
        youngs : float, optional (Pa)
            Young's modulus.
        method : {'airy', 'flexure'}, optional
            Method to use to calculate deflections.
        rho_mantle : float, optional
            Density of the mantle (kg / m^3).
        rho_water : float, optional
            Density of the sea water (kg / m^3).
        gravity : float, optional
            Acceleration due to gravity (m / s^2).
        rows : int, optional
            Node rows that this component will operate on (default is to
            operate on *all* rows).
        """
        if method not in ("airy", "flexure"):
            raise ValueError(f"{method}: method not understood")

        super().__init__(grid)

        self._method = method
        self.youngs = youngs
        self.rho_mantle = rho_mantle
        self.rho_water = rho_water
        self.gravity = gravity
        self.eet = eet

        self.initialize_output_fields()

        self._rows = (rows,) or Ellipsis

        self._x_at_node = self._grid.x_of_node.reshape(self._grid.shape).copy()

    @property
    def eet(self):
        """Effective elastic thickness (m)."""
        return self._eet

    @eet.setter
    def eet(self, new_val):
        if new_val <= 0:
            raise ValueError("Effective elastic thickness must be positive.")
        with contextlib.suppress(AttributeError):
            del self._rigidity, self._alpha
        self._eet = float(new_val)

    @property
    def youngs(self):
        """Young's modulus of lithosphere (Pa)."""
        return self._youngs

    @youngs.setter
    def youngs(self, new_val):
        if new_val <= 0:
            raise ValueError("Young's modulus must be positive.")
        with contextlib.suppress(AttributeError):
            del self._rigidity, self._alpha
        self._youngs = float(new_val)

    @property
    def rho_water(self):
        """Density of water (kg/m^3)."""
        return self._rho_water

    @rho_water.setter
    def rho_water(self, new_val):
        if new_val <= 0:
            raise ValueError("Water density must be positive.")
        with contextlib.suppress(AttributeError):
            del self._gamma_mantle, self._alpha
        self._rho_water = float(new_val)

    @property
    def rho_mantle(self):
        """Density of mantle (kg/m^3)."""
        return self._rho_mantle

    @rho_mantle.setter
    def rho_mantle(self, new_val):
        if new_val <= 0:
            raise ValueError("Mantle density must be positive.")
        with contextlib.suppress(AttributeError):
            del self._gamma_mantle, self._alpha
        self._rho_mantle = float(new_val)

    @property
    def gravity(self):
        """Acceleration due to gravity (m/s^2)."""
        return self._gravity

    @gravity.setter
    def gravity(self, new_val):
        if new_val <= 0:
            raise ValueError("Acceleration due to gravity must be positive.")
        with contextlib.suppress(AttributeError):
            del self._gamma_mantle, self._alpha
        self._gravity = float(new_val)

    @property
    def alpha(self):
        """Flexure parameter (m)."""
        try:
            self._alpha
        except AttributeError:
            self._alpha = np.power(4 * self.rigidity / self.gamma_mantle, 0.25)
        return self._alpha

    @property
    def rigidity(self):
        """Flexural rigidity (N m)."""
        try:
            self._rigidity
        except AttributeError:
            self._rigidity = (
                self._eet**3.0 * self._youngs / (12.0 * (1.0 - self._POISSON**2.0))
            )
        return self._rigidity

    @property
    def gamma_mantle(self):
        """Specific density of mantle (N/m^3)."""
        try:
            self._gamma_mantle
        except AttributeError:
            self._gamma_mantle = (self._rho_mantle - self._rho_water) * self._gravity
        return self._gamma_mantle

    @property
    def method(self):
        """Name of method used to calculate deflections."""
        return self._method

    @property
    def x_at_node(self):
        return self._x_at_node

    @property
    def load_at_node(self):
        return self._grid.at_node[
            "lithosphere__increment_of_overlying_pressure"
        ].reshape(self._grid.shape)

    @property
    def dz_at_node(self):
        return self._grid.at_node[
            "lithosphere_surface__increment_of_elevation"
        ].reshape(self._grid.shape)

    def update(self):
        """Update fields with current loading conditions."""
        load = self.load_at_node[self._rows]
        deflection = self.dz_at_node[self._rows]

        if self._method == "airy":
            deflection[:] = load / self.gamma_mantle
        else:
            Flexure1D.calc_flexure(
                self.x_at_node[0], load, self.alpha, self.rigidity, out=deflection
            )

        if not np.may_share_memory(deflection, self.dz_at_node):
            self.dz_at_node[self._rows] = deflection

    def run_one_step(self):
        self.update()

    def subside_loads(self, loads, out=None):
        """Subside surface due to multiple loads.

        Parameters
        ----------
        loads : ndarray of float
            Loads applied to each grid node (Pa).
        out : ndarray of float, optional
            Buffer to place resulting deflection values (m).

        Returns
        -------
        ndarray of float
            Deflections caused by the loading.
        """
        if out is None:
            out = np.zeros(self._grid.shape)
        loads = np.asarray(loads)
        if self._method == "airy":
            out[:] = loads / self.gamma_mantle
        else:
            Flexure1D.calc_flexure(
                self.x_at_node[0], loads, self.alpha, self.rigidity, out=out
            )

        return out

    @staticmethod
    def calc_flexure(x, loads, alpha, rigidity, out=None):
        """Subside surface due to multiple loads.

        Parameters
        ----------
        x : ndarray of float
            Position of grid nodes (m).
        loads : ndarray of float
            Loads applied to each grid node (Pa).
        alpha : float
            Flexure parameter of the lithosphere (m).
        rigidity : float
            Flexural rigidity of the lithosphere (N m).
        out : ndarray of float, optional
            Buffer to place resulting deflection values (m).

        Returns
        -------
        ndarray of float
            Deflections caused by the loading.
        """
        if out is None:
            out = np.zeros_like(loads, dtype=float)

        loads = loads.reshape((-1, loads.shape[-1]))
        dz = out[..., :].reshape((-1, out.shape[-1]))
        subside_load_1d(x, loads, alpha, rigidity, dz[..., :])

        return out

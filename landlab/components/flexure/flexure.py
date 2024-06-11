#!/usr/bin/env python
"""Deform the lithosphere with 1D or 2D flexure.

Landlab component that implements a 1 and 2D lithospheric flexure
model.

Examples
--------

Create a grid on which we will run the flexure calculations.

>>> from landlab import RasterModelGrid
>>> from landlab.components.flexure import Flexure
>>> grid = RasterModelGrid((5, 4), xy_spacing=(1.0e4, 1.0e4))
>>> lith_press = grid.add_zeros("lithosphere__overlying_pressure_increment", at="node")

Check the fields that are used as input to the flexure component.

>>> Flexure.input_var_names
('lithosphere__overlying_pressure_increment',)

Check the units for the fields.

>>> Flexure.var_units("lithosphere__overlying_pressure_increment")
'Pa'

If you are not sure about one of the input or output variables, you can
get help for specific variables.

>>> Flexure.var_help("lithosphere__overlying_pressure_increment")
name: lithosphere__overlying_pressure_increment
description:
  Applied pressure to the lithosphere over a time step
units: Pa
unit agnostic: True
at: node
intent: in

>>> flex = Flexure(grid)

In creating the component, a field (initialized with zeros) was added to the
grid. Reset the interior nodes for the loading.

>>> dh = grid.at_node["lithosphere__overlying_pressure_increment"]
>>> dh = dh.reshape(grid.shape)
>>> dh[1:-1, 1:-1] = flex.gamma_mantle

>>> flex.update()

>>> flex.output_var_names
('lithosphere_surface__elevation_increment',)
>>> flex.grid.at_node["lithosphere_surface__elevation_increment"].reshape(grid.shape)
array([[0., 0., 0., 0.],
       [0., 1., 1., 0.],
       [0., 1., 1., 0.],
       [0., 1., 1., 0.],
       [0., 0., 0., 0.]])
"""

import numpy as np

from landlab import Component
from landlab.components.flexure._ext.flexure2d import subside_loads as _subside_loads
from landlab.components.flexure._ext.flexure2d_slow import subside_grid_in_parallel
from landlab.components.flexure.funcs import get_flexure_parameter


class Flexure(Component):
    """Deform the lithosphere with 1D or 2D flexure.

    Landlab component that implements a 1 and 2D lithospheric flexure
    model.

    Examples
    --------

    >>> from landlab import RasterModelGrid
    >>> from landlab.components.flexure import Flexure
    >>> grid = RasterModelGrid((5, 4), xy_spacing=(1.0e4, 1.0e4))
    >>> lith_press = grid.add_zeros(
    ...     "lithosphere__overlying_pressure_increment", at="node"
    ... )

    >>> flex = Flexure(grid)
    >>> flex.name
    'Flexure'
    >>> flex.input_var_names
    ('lithosphere__overlying_pressure_increment',)
    >>> flex.output_var_names
    ('lithosphere_surface__elevation_increment',)
    >>> sorted(flex.units)
    [('lithosphere__overlying_pressure_increment', 'Pa'),
     ('lithosphere_surface__elevation_increment', 'm')]

    >>> flex.grid.number_of_node_rows
    5
    >>> flex.grid.number_of_node_columns
    4
    >>> flex.grid is grid
    True

    >>> np.all(grid.at_node["lithosphere_surface__elevation_increment"] == 0.0)
    True

    >>> np.all(grid.at_node["lithosphere__overlying_pressure_increment"] == 0.0)
    True
    >>> flex.update()
    >>> np.all(grid.at_node["lithosphere_surface__elevation_increment"] == 0.0)
    True

    >>> load = grid.at_node["lithosphere__overlying_pressure_increment"]
    >>> load[4] = 1e9
    >>> dz = grid.at_node["lithosphere_surface__elevation_increment"]
    >>> np.all(dz == 0.0)
    True

    >>> flex.update()
    >>> np.all(grid.at_node["lithosphere_surface__elevation_increment"] == 0.0)
    False

    References
    ----------
    **Required Software Citation(s) Specific to this Component**

    Hutton, E., Syvitski, J. (2008). Sedflux 2.0: An advanced process-response
    model that generates three-dimensional stratigraphy. Computers &
    Geosciences.  34(10), 1319-1337.
    https://dx.doi.org/10.1016/j.cageo.2008.02.013

    **Additional References**

    Lambeck, K.: Geophysical Geodesy, The Slow Deformations of the Earth,
    Clarendon Press, Oxford, UK, 718 pp., 1988.

    """

    _name = "Flexure"

    _unit_agnostic = True

    _cite_as = r"""
    @article{hutton2008sedflux,
        title={Sedflux 2.0: An advanced process-response model that generates
               three-dimensional stratigraphy},
        author={Hutton, Eric WH and Syvitski, James PM},
        journal={Computers \& Geosciences},
        volume={34},
        number={10},
        pages={1319--1337},
        year={2008},
        publisher={Pergamon}
        }"""

    _info = {
        "lithosphere__overlying_pressure_increment": {
            "dtype": float,
            "intent": "in",
            "optional": False,
            "units": "Pa",
            "mapping": "node",
            "doc": "Applied pressure to the lithosphere over a time step",
        },
        "lithosphere_surface__elevation_increment": {
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

    def __init__(
        self,
        grid,
        eet=65e3,
        youngs=7e10,
        method="airy",
        rho_mantle=3300.0,
        gravity=9.80665,
        n_procs=1,
    ):
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
        n_procs : int, optional
            Number of processors to use for calculations.
        """
        if method not in ("airy", "flexure"):
            raise ValueError(f"{method}: method not understood")

        super().__init__(grid)

        self._youngs = youngs
        self._method = method
        self._rho_mantle = rho_mantle
        self._gravity = gravity
        self.eet = eet
        self._n_procs = n_procs

        self.initialize_output_fields()

        self._r = self._create_kei_func_grid(
            self._grid.shape, (self._grid.dy, self._grid.dx), self.alpha
        )

    @property
    def eet(self):
        """Effective elastic thickness (m)."""
        return self._eet

    @eet.setter
    def eet(self, new_val):
        if new_val <= 0:
            raise ValueError("Effective elastic thickness must be positive.")
        self._eet = new_val
        self._r = self._create_kei_func_grid(
            self._grid.shape, (self._grid.dy, self._grid.dx), self.alpha
        )

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
        return get_flexure_parameter(
            self._eet, self._youngs, 2, gamma_mantle=self.gamma_mantle
        )

    @staticmethod
    def _create_kei_func_grid(shape, xy_spacing, alpha):
        from scipy.special import kei

        dx, dy = np.meshgrid(
            np.arange(shape[1]) * xy_spacing[1], np.arange(shape[0]) * xy_spacing[0]
        )

        return kei(np.sqrt(dx**2 + dy**2) / alpha)

    def update(self):
        """Update fields with current loading conditions."""
        load = self._grid.at_node["lithosphere__overlying_pressure_increment"]
        deflection = self._grid.at_node["lithosphere_surface__elevation_increment"]

        new_load = load.copy()

        deflection.fill(0.0)

        if self.method == "airy":
            deflection[:] = new_load / self.gamma_mantle
        else:
            self.subside_loads(new_load, out=deflection)

    def subside_loads_slow(self, loads, out=None):
        """Subside surface due to multiple loads.

        Parameters
        ----------
        loads : ndarray of float
            Loads applied to each grid node.
        out : ndarray of float, optional
            Buffer to place resulting deflection values.

        Returns
        -------
        ndarray of float
            Deflections caused by the loading.
        """
        if out is None:
            out = np.zeros(self._grid.shape, dtype=float)
        dz = out.reshape(self._grid.shape)
        load = loads.reshape(self._grid.shape)

        subside_grid_in_parallel(
            dz,
            load * self._grid.dx * self._grid.dy,
            self._r,
            self.alpha,
            self.gamma_mantle,
            self._n_procs,
        )

        return out

    def subside_loads(self, loads, row_col_of_load=None, out=None):
        """Subside surface due to multiple loads.

        Parameters
        ----------
        loads : ndarray of float
            Loads applied to grid node. ``loads`` can be either an array
            of size ``n_nodes``, in which case the load values are applied
            at their corresponding nodes, or an array of arbitray length,
            in which case loads are applied at locations supplied through
            the ``row_col_of_load`` keyword.
        row_col_of_load: tuple of ndarray of int, optional
            If provided, the row and column indices where loads are applied.
            The first element of the tuple is an array of rows while the
            seconds element is an array of columns.
        out : ndarray of float, optional
            Buffer to place resulting deflection values. If not provided,
            deflections will be placed into a newly-allocated array.

        Returns
        -------
        ndarray of float
            Deflections caused by the loading.
        """
        loads = np.asarray(loads)
        if out is None:
            out = self.grid.zeros(at="node")
        dz = out.reshape(self.grid.shape)

        if row_col_of_load is None:
            loads, row_col_of_load = self._handle_no_row_col(loads)
        row_of_load, col_of_load = row_col_of_load

        _subside_loads(
            dz,
            self._r.reshape(self.grid.shape),
            loads * self.grid.dx * self.grid.dy,
            np.asarray(row_of_load),
            np.asarray(col_of_load),
            self.alpha,
            self.gamma_mantle,
        )

        return out

    def _handle_no_row_col(self, loads, tol=1e-6):
        """Handle the case where the row_col_of_load keyword is not provided."""
        loads = loads.reshape(self.grid.shape)
        row_col_of_load = np.unravel_index(
            np.flatnonzero(np.abs(loads) > tol), self.grid.shape
        )
        loads = loads[row_col_of_load]

        return loads, row_col_of_load

# -*- coding: utf-8 -*-
"""Landlab component for overland flow using the diffusion-wave approximation.

Created on Fri May 27 14:26:13 2016

@author: gtucker
"""


import numpy as np

from landlab import Component

_FOUR_THIRDS = 4.0 / 3.0


class LinearDiffusionOverlandFlowRouter(Component):
    """Calculate water flow over topography.

    Landlab component that implements a two-dimensional, linearized
    diffusion-wave model. The diffusion-wave approximation is a simplification
    of the shallow-water equations that omits the momentum terms. The flow
    velocity is calculated using the local water-surface slope as an
    approximation of the energy slope. With this linearized form, flow velocity
    is calculated using a linearized Manning equation, with the water-surface
    slope being used as the slope factor. There are two governing equations, one
    that represents conservation of water mass:

    ..math::

        \frac{\partial H}{\partial t} = (P - I) - \nabla\cdot\mathbf{q}

    where :math:`H(x,y,t)` is local water depth, :math:`t` is time, :math:`P`
    is precipitation rate, :math:`I` is infiltration rate, and :math:`\mathbf{q}`
    is specific water discharge, which equals depth times depth-averaged
    velocity. The other governing equation represents specific discharge as a
    function of gravity, pressure, and friction:

    ..math::

        \mathbf{q} = \frac{H^{4/3}}{n^2 U_c} \nabla w

    where :math:`n` is the friction factor ("Manning's n"), :math:`U_c` is a
    characteristic flow velocity, and :math:`w` is water-surface height, which
    is the sum of topographic elevation plus water depth.

    Examples
    --------
    >>> from landlab import RasterModelGrid
    >>> grid = RasterModelGrid((3, 3))
    >>> _ = grid.add_zeros('topographic__elevation', at='node')
    >>> olf = LinearDiffusionOverlandFlowRouter(grid, roughness=0.1)
    >>> round(olf.vel_coef)
    100
    >>> olf.precip_rate
    1e-05
    >>> olf.precip_rate = 1.0e-4
    >>> olf.run_one_step(dt=10.0)
    >>> grid.at_node['surface_water__depth'][4]
    0.001

    References
    ----------
    **Required Software Citation(s) Specific to this Component**

    None Listed

    **Additional References**

    None Listed

    """

    _name = "LinearDiffusionOverlandFlowRouter"

    _unit_agnostic = True

    _info = {
        "surface_water__depth": {
            "dtype": float,
            "intent": "out",
            "optional": False,
            "units": "m",
            "mapping": "node",
            "doc": "Depth of water on the surface",
        },
        "topographic__elevation": {
            "dtype": float,
            "intent": "in",
            "optional": False,
            "units": "m",
            "mapping": "node",
            "doc": "Land surface topographic elevation",
        },
        "water__specific_discharge": {
            "dtype": float,
            "intent": "out",
            "optional": False,
            "units": "m2/s",
            "mapping": "link",
            "doc": "flow discharge component in the direction of the link",
        },
        "water__velocity": {
            "dtype": float,
            "intent": "out",
            "optional": False,
            "units": "m/s",
            "mapping": "link",
            "doc": "flow velocity component in the direction of the link",
        },
        "water_surface__gradient": {
            "dtype": float,
            "intent": "out",
            "optional": False,
            "units": "m/s",
            "mapping": "link",
            "doc": "gradient of water surface",
        },
    }

    def __init__(
        self,
        grid,
        precip_rate=1.0e-5,
        infilt_rate=0.0,
        roughness=0.01,
        velocity_scale=1.0,
    ):
        """Initialize the LinearDiffusionOverlandFlowRouter.

        Parameters
        ----------
        grid : ModelGrid
            Landlab ModelGrid object
        precip_rate : float, optional (defaults to 36 mm/hr)
            Precipitation rate, m/s
        infilt_rate : float, optional (defaults to 0)
            Maximum rate of infiltration, m/s
        roughness : float, defaults to 0.01
            Manning roughness coefficient, s/m^1/3
        velocity_scale : float, defaults to 1
            Characteristic flow velocity, m/s
        """
        super().__init__(grid)

        # Store parameters and do unit conversion
        self._precip = precip_rate
        self._infilt = infilt_rate
        self._vel_coef = 1.0 / (roughness ** 2 * velocity_scale)

        self._elev = grid.at_node["topographic__elevation"]

        self.initialize_output_fields()
        self._depth = grid.at_node["surface_water__depth"]
        self._vel = grid.at_link["water__velocity"]
        self._disch = grid.at_link["water__specific_discharge"]
        self._wsgrad = grid.at_link["water_surface__gradient"]

    @property
    def precip_rate(self):
        """Precipitation rate"""
        return self._precip

    @precip_rate.setter
    def precip_rate(self, value):
        self._precip = value

    @property
    def vel_coef(self):
        """Velocity coefficient.

        (1/(roughness^2 x velocity_scale)
        """
        return self._vel_coef

    def run_one_step(self, dt):
        """Calculate water flow for a time period `dt`.

        Default units for dt are *seconds*.
        """
        # Calculate water depth at links. This implements an "upwind" scheme
        # in which water depth at the links is the depth at the higher of the
        # two nodes.
        H_link = self._grid.map_value_at_max_node_to_link(
            "topographic__elevation", "surface_water__depth"
        )

        # Calculate the water-surface gradient
        self.grid.calc_grad_at_link(self._elev + self._depth, out=self._wsgrad)

        # Calculate velocity using the linearized Manning equation.
        self._vel = -self._vel_coef * H_link ** _FOUR_THIRDS * self._wsgrad

        # Calculate discharge
        self._disch[:] = H_link * self._vel

        # Flux divergence
        dqda = self._grid.calc_flux_div_at_node(self._disch)

        # Rates of infiltration and runoff (TODO: do better than this hack)
        infilt = np.minimum(self._infilt, self._depth / dt)

        # Rate of change of water depth
        dHdt = self._precip - infilt - dqda

        # Update water depth: simple forward Euler scheme
        self._depth[self._grid.core_nodes] += dHdt[self._grid.core_nodes] * dt

        # Very crude numerical hack: prevent negative water depth (TODO: better)
        self._depth[np.where(self._depth < 0.0)[0]] = 0.0

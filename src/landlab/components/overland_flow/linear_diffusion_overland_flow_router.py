"""Landlab component for overland flow using the linearized diffusion-wave approximation.

Created on Fri May 27 14:26:13 2016

@author: gtucker
"""

import numpy as np

from landlab import Component

_FOUR_THIRDS = 4.0 / 3.0
_SEVEN_THIRDS = 7.0 / 3.0
_MICRO_DEPTH = 1.0e-6  # tiny water depth to avoid blowup in time-step estimator


class LinearDiffusionOverlandFlowRouter(Component):
    r"""Calculate water flow over topography.

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

    Infiltration rate should decline smoothly to zero as surface water depth
    approaches zero. To ensure this, infiltration rate is calculated as

    ..math::

        I = I_c \left( 1 - e^{-H / H_i} ) \right)

    where :math:`H_i` is a characteristic water depth. The concept here is that
    when :math:`H \le H_i`, small spatial variations in water depth will leave
    parts of the ground un-ponded and therefore not subject to any infiltration.
    Mathematically, :math:`I \approx 0.95 I_c` when :math:`H = 3H_i`.

    Examples
    --------
    >>> from landlab import RasterModelGrid
    >>> grid = RasterModelGrid((3, 3))
    >>> _ = grid.add_zeros("topographic__elevation", at="node")
    >>> olf = LinearDiffusionOverlandFlowRouter(grid, roughness=0.1)
    >>> round(olf.vel_coef)
    100
    >>> olf.rain_rate
    1e-05
    >>> olf.rain_rate = 1.0e-4
    >>> olf.run_one_step(dt=10.0)
    >>> grid.at_node["surface_water__depth"][4]
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
        "surface_water__depth_at_link": {
            "dtype": float,
            "intent": "out",
            "optional": False,
            "units": "m",
            "mapping": "link",
            "doc": "Depth of water on the surface at grid links",
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
            "doc": "Downstream gradient of the water surface.",
        },
        "water_surface__elevation": {
            "dtype": float,
            "intent": "out",
            "optional": False,
            "units": "m",
            "mapping": "node",
            "doc": "Elevation of the water surface.",
        },
    }

    def __init__(
        self,
        grid,
        roughness=0.01,
        rain_rate=1.0e-5,
        infilt_rate=0.0,
        infilt_depth_scale=0.001,
        velocity_scale=1.0,
        cfl_factor=0.2,
    ):
        """Initialize the LinearDiffusionOverlandFlowRouter.

        Parameters
        ----------
        grid : ModelGrid
            Landlab ModelGrid object
        roughness : float, defaults to 0.01
            Manning roughness coefficient, s/m^1/3
        rain_rate : float, optional (defaults to 36 mm/hr)
            Rainfall rate, m/s
        infilt_depth_scale : float, defaults to 0.001
            Depth scale for water infiltration, m
        infilt_rate : float, optional (defaults to 0)
            Maximum rate of infiltration, m/s
        velocity_scale : float, defaults to 1
            Characteristic flow velocity, m/s
        cfl_factor : float, optional (defaults to 0.2)
            Time-step control factor: fraction of maximum estimated time-step
            that is actually used (must be <=1)
        """
        super().__init__(grid)

        if roughness <= 0.0:
            raise ValueError(f"roughness must be greater than zero {roughness}")
        if rain_rate < 0.0:
            raise ValueError(f"rain_rate must not be negative {rain_rate}")
        if infilt_depth_scale <= 0.0:
            raise ValueError(
                f"infilt_depth_scale must be greater than zero {infilt_depth_scale}"
            )
        if infilt_rate < 0.0:
            raise ValueError(f"infilt_rate must not be negative {infilt_rate}")
        if velocity_scale <= 0.0:
            raise ValueError(
                f"velocity_scale must be greater than zero {velocity_scale}"
            )
        if cfl_factor > 1.0 or cfl_factor <= 0.0:
            raise ValueError(f"cfl_factor must >0, <=1 {cfl_factor}")

        # Store parameters and do unit conversion
        self._rain = rain_rate
        self._infilt = infilt_rate
        self._infilt_depth_scale = infilt_depth_scale
        self._vel_coef = 1.0 / (roughness**2 * velocity_scale)

        self._elev = grid.at_node["topographic__elevation"]

        self.initialize_output_fields()
        self._depth = grid.at_node["surface_water__depth"]
        self._depth_at_link = grid.at_link["surface_water__depth_at_link"]
        self._vel = grid.at_link["water__velocity"]
        self._disch = grid.at_link["water__specific_discharge"]
        self._wsgrad = grid.at_link["water_surface__gradient"]
        self._water_surf_elev = grid.at_node["water_surface__elevation"]

        self._inactive_links = grid.status_at_link == grid.BC_LINK_IS_INACTIVE

        self._cfl_param = cfl_factor * 0.5 * np.amin(grid.length_of_link) ** 2

    @property
    def rain_rate(self):
        """Rainfall rate"""
        return self._rain

    @rain_rate.setter
    def rain_rate(self, value):
        if value < 0.0:
            raise ValueError(f"rain_rate must be positive {value}")
        self._rain = value

    @property
    def vel_coef(self):
        """Velocity coefficient.

        (1/(roughness^2 x velocity_scale)
        """
        return self._vel_coef

    def _cfl_time_step(self):
        """Calculate maximum time-step size using CFL criterion for explicit
        FTCS diffusion."""
        max_water_depth = np.amax(self._depth, initial=_MICRO_DEPTH)
        max_diffusivity = self._vel_coef * max_water_depth**_SEVEN_THIRDS
        return self._cfl_param / max_diffusivity

    def update_for_one_iteration(self, iter_dt):
        """Update state variables for one iteration of duration iter_dt."""

        # Calculate the water-surface elevation
        self._water_surf_elev[:] = self._elev + self._depth

        # Calculate water depth at links. This implements an "upwind" scheme
        # in which water depth at the links is the depth at the higher of the
        # two nodes.
        self._grid.map_value_at_max_node_to_link(
            self._water_surf_elev, "surface_water__depth", out=self._depth_at_link
        )

        # Calculate the water-surface gradient and impose any closed boundaries
        self.grid.calc_grad_at_link(self._water_surf_elev, out=self._wsgrad)
        self._wsgrad[self._inactive_links] = 0.0

        # Calculate velocity using the linearized Manning equation.
        self._vel[:] = (
            -self._vel_coef * self._depth_at_link**_FOUR_THIRDS * self._wsgrad
        )

        # Calculate discharge
        self._disch[:] = self._depth_at_link * self._vel

        # Flux divergence
        dqda = self._grid.calc_flux_div_at_node(self._disch)

        # Rates of infiltration and runoff
        infilt = self._infilt * (1.0 - np.exp(-self._depth / self._infilt_depth_scale))

        # Rate of change of water depth
        dHdt = self._rain - infilt - dqda

        # Update water depth: simple forward Euler scheme
        self._depth[self._grid.core_nodes] += dHdt[self._grid.core_nodes] * iter_dt

        # Very crude numerical hack: prevent negative water depth (TODO: better)
        self._depth.clip(min=0.0, out=self._depth)

    def run_one_step(self, dt):
        """Calculate water flow for a time period `dt`.

        Default units for dt are *seconds*. We use a time-step subdivision
        algorithm that ensures step size is always below CFL limit.
        """
        remaining_time = dt
        while remaining_time > 0.0:
            dtmax = self._cfl_time_step()  # biggest possible time step size
            dt_this_iter = min(dtmax, remaining_time)  # step we'll actually use
            self.update_for_one_iteration(dt_this_iter)
            remaining_time -= dt_this_iter

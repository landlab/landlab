"""Landlab component for overland flow using the kinematic-wave approximation.

Created on Fri May 27 14:26:13 2016

@author: gtucker
"""

import numpy as np

from landlab import Component


class KinwaveOverlandFlowModel(Component):
    """Calculate water flow over topography.

    Landlab component that implements a two-dimensional
    kinematic wave model. This is an extremely simple, unsophisticated
    model, originally built simply to demonstrate the component creation
    process. Limitations to the present version include: infiltration is
    handled very crudely, the called is responsible for picking a stable
    time step size (no adaptive time stepping is used in the `run_one_step`
    method), precipitation rate is constant for a given duration (then zero),
    and all parameters are uniform in space. Also, the terrain is assumed
    to be stable over time. Caveat emptor!

    Examples
    --------
    >>> from landlab import RasterModelGrid
    >>> rg = RasterModelGrid((4, 5), xy_spacing=10.0)
    >>> z = rg.add_zeros("topographic__elevation", at="node")
    >>> s = rg.add_zeros("topographic__gradient", at="link")
    >>> kw = KinwaveOverlandFlowModel(rg)
    >>> kw.vel_coef
    100.0
    >>> rg.at_node["surface_water__depth"]
    array([0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.,
           0., 0., 0., 0., 0., 0., 0.])

    References
    ----------
    **Required Software Citation(s) Specific to this Component**

    None Listed

    **Additional References**

    None Listed

    """

    _name = "KinwaveOverlandFlowModel"

    _unit_agnostic = False

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
        "topographic__gradient": {
            "dtype": float,
            "intent": "in",
            "optional": False,
            "units": "m/m",
            "mapping": "link",
            "doc": "Gradient of the ground surface",
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
    }

    def __init__(
        self,
        grid,
        precip_rate=1.0,
        precip_duration=1.0,
        infilt_rate=0.0,
        roughness=0.01,
    ):
        """Initialize the KinwaveOverlandFlowModel.

        Parameters
        ----------
        grid : ModelGrid
            Landlab ModelGrid object
        precip_rate : float, optional (defaults to 1 mm/hr)
            Precipitation rate, mm/hr
        precip_duration : float, optional (defaults to 1 hour)
            Duration of precipitation, hours
        infilt_rate : float, optional (defaults to 0)
            Maximum rate of infiltration, mm/hr
        roughness : float, defaults to 0.01
            Manning roughness coefficient, s/m^1/3
        """
        super().__init__(grid)

        # Store parameters and do unit conversion
        self._current_time = 0

        self._precip = precip_rate / 3600000.0  # convert to m/s
        self._precip_duration = precip_duration * 3600.0  # h->s
        self._infilt = infilt_rate / 3600000.0  # convert to m/s
        self._vel_coef = 1.0 / roughness  # do division now to save time

        # Create fields...
        #   Elevation
        self._elev = grid.at_node["topographic__elevation"]

        #   Slope
        self._slope = grid.at_link["topographic__gradient"]

        self.initialize_output_fields()
        self._depth = grid.at_node["surface_water__depth"]
        self._vel = grid.at_link["water__velocity"]
        self._disch = grid.at_link["water__specific_discharge"]

        # Calculate the ground-surface slope (assume it won't change)
        self._slope[self._grid.active_links] = self._grid.calc_grad_at_link(self._elev)[
            self._grid.active_links
        ]
        self._sqrt_slope = np.sqrt(self._slope)
        self._sign_slope = np.sign(self._slope)

    @property
    def vel_coef(self):
        """Velocity coefficient.

        (1/roughness)
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

        # Calculate velocity using the Manning equation.
        self._vel = (
            -self._sign_slope * self._vel_coef * H_link**0.66667 * self._sqrt_slope
        )

        # Calculate discharge
        self._disch[:] = H_link * self._vel

        # Flux divergence
        dqda = self._grid.calc_flux_div_at_node(self._disch)

        # Rate of change of water depth
        if self._current_time < self._precip_duration:
            ppt = self._precip
        else:
            ppt = 0.0
        dHdt = ppt - self._infilt - dqda

        # Update water depth: simple forward Euler scheme
        self._depth[self._grid.core_nodes] += dHdt[self._grid.core_nodes] * dt

        # Very crude numerical hack: prevent negative water depth
        self._depth[np.where(self._depth < 0.0)[0]] = 0.0

        self._current_time += dt


if __name__ == "__main__":
    import doctest

    doctest.testmod()

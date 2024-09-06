#! /usr/bin/env python
from contextlib import contextmanager

import numpy as np

from landlab import Component

_EPSILON = 1.0e-6  # tiny negative depth avoid fast carb growth at zero depth


@contextmanager
def set_numpy_err(*args, **kwds):
    settings = np.seterr(*args, **kwds)
    try:
        yield settings
    finally:
        np.seterr(**settings)


def smooth_heaviside(x, width=0.5, out=None):
    """Return a smoothed heaviside function (step function).

    Parameters
    ----------
    x : array or float
        Dependent variable
    width : float (default 0.5)
        Width parameter for smoothing (same units as x)
    out : array (default None)
        Optional array in which to store result; must have len(x)

    Examples
    --------
    >>> import numpy as np
    >>> np.round(smooth_heaviside(np.array([-1, 0, 1])), 3)
    array([0.018, 0.5  , 0.982])
    >>> smooth_heaviside(np.array([-1, 0, 1]), width=0.0)
    array([0. , 0.5, 1. ])
    """
    if width > 0.0:
        with set_numpy_err(over="ignore"):
            return np.divide(1.0, (1.0 + np.exp(-2.0 / width * np.asarray(x))), out=out)
    else:
        return np.heaviside(x, 0.5)


class CarbonateProducer(Component):
    """Calculate marine carbonate production and deposition.

    Uses the growth-rate of Bosscher and Schlager (1992), which represents the
    vertical growth rate :math:`G` as:

    ..math::

        G = G_m \tanh (I_0 e^{-kz} / I_k)

    where :math:`G_m` is the maximum growth rate, :math:`I_0` is the surface
    light intensity, :math:`I_k` is the saturating light intensity, :math:`z` is
    water depth, and :math:`k` is the light extinction coefficient.

    Bosscher and Schlager (1992) suggest plausible values or ranges for these
    parameters as follows: :math:`G_m` = 10 to 15 mm/y, :math:`I_0` = 2,000 to
    2,250 micro Einsteins per square meter per second in the tropics, :math:`I_k`
    = 50 to 450 micro Einsteins per square meter per second, and :math:`k` =
    0.04 to 0.16 m:math:`^{-1}` (corresponding to a decay depth of 6.25 to 16 m).
    The default values used in this component are based on these estimates.

    Examples
    --------
    >>> from landlab import RasterModelGrid
    >>> from landlab.components import CarbonateProducer
    >>> grid = RasterModelGrid((3, 3))
    >>> elev = grid.add_zeros("topographic__elevation", at="node")
    >>> sealevel = grid.add_field("sea_level__elevation", 0.0, at="grid")
    >>> elev[:] = -1.0
    >>> cp = CarbonateProducer(grid)
    >>> np.round(cp.calc_carbonate_production_rate()[4], 2)
    0.01
    >>> elev[:] = -40.0
    >>> np.round(cp.calc_carbonate_production_rate()[4], 5)
    0.00091
    >>> cp.sea_level = -39.0
    >>> np.round(cp.calc_carbonate_production_rate()[4], 2)
    0.01
    >>> thickness = cp.produce_carbonate(10.0)
    >>> np.round(10 * grid.at_node["carbonate_thickness"][4])
    1.0
    >>> thickness is grid.at_node["carbonate_thickness"]
    True
    >>> cp.run_one_step(10.0)
    >>> np.round(10 * thickness[4])
    2.0

    References
    ----------
    Bosscher, H., & Schlager, W. (1992). Computer simulation of reef growth.
    Sedimentology, 39(3), 503-512.
    """

    _name = "CarbonateProducer"

    _unit_agnostic = True

    _info = {
        "sea_level__elevation": {
            "dtype": "float",
            "intent": "in",
            "optional": False,
            "units": "m",
            "mapping": "grid",
            "doc": "Sea level elevation",
        },
        "topographic__elevation": {
            "dtype": "float",
            "intent": "inout",
            "optional": False,
            "units": "m",
            "mapping": "node",
            "doc": "Land surface topographic elevation",  # and seafloor
        },
        "carbonate_production_rate": {
            "dtype": "float",
            "intent": "out",
            "optional": False,
            "units": "m / y",
            "mapping": "node",
            "doc": "Carbonate production rate",
        },
        "carbonate_thickness": {
            "dtype": "float",
            "intent": "out",
            "optional": True,
            "units": "m",
            "mapping": "node",
            "doc": "Carbonate thickness",
        },
        "water_depth": {
            "dtype": "float",
            "intent": "out",
            "optional": False,
            "units": "m",
            "mapping": "node",
            "doc": "Water depth",
        },
    }

    def __init__(
        self,
        grid,
        max_carbonate_production_rate=0.01,
        extinction_coefficient=0.1,
        surface_light=2000.0,
        saturating_light=400.0,
        tidal_range=0.0,
    ):
        """
        Parameters
        ----------
        grid: ModelGrid (RasterModelGrid, HexModelGrid, etc.)
            A landlab grid object.
        max_carbonate_production_rate: float (default 0.01 m/y)
            Maximum rate of carbonate production, in m thickness per year
        extinction_coefficient: float (default 0.1 m^-1)
            Coefficient of light extinction in water column
        surface_light: float (default 2000.0 micro Einstein per m2 per s)
            Light intensity (photosynthetic photon flux density) at sea surface.
        saturating_light: float (default 400.0 micro Einstein per m2 per s)
            Saturating light intensity (photosynthetic photon flux density)
        tidal_range: float (default zero)
            Tidal range used to smooth growth rate when surface is near mean
            sea level.
        """
        super().__init__(grid)

        self.extinction_coefficient = extinction_coefficient
        self.max_carbonate_production_rate = max_carbonate_production_rate
        self.surface_light = surface_light
        self.tidal_range = tidal_range
        self.saturating_light = saturating_light

        super().initialize_output_fields()
        self._carb_prod_rate = grid.at_node["carbonate_production_rate"]
        self._depth = grid.at_node["water_depth"]

        if "carbonate_thickness" not in grid.at_node:
            grid.add_zeros("carbonate_thickness", at="node")
        self._carbonate_thickness = grid.at_node["carbonate_thickness"]

    @property
    def extinction_coefficient(self):
        return self._extinction_coefficient

    @extinction_coefficient.setter
    def extinction_coefficient(self, value):
        if value > 0.0:
            self._extinction_coefficient = float(value)
        else:
            raise ValueError("must be non-negative")

    @property
    def max_carbonate_production_rate(self):
        return self._max_carbonate_production_rate

    @max_carbonate_production_rate.setter
    def max_carbonate_production_rate(self, value):
        if value >= 0.0:
            self._max_carbonate_production_rate = float(value)
        else:
            raise ValueError("must be non-negative")

    @property
    def surface_light(self):
        return self._surface_light

    @surface_light.setter
    def surface_light(self, value):
        if value >= 0.0:
            self._surface_light = float(value)
        else:
            raise ValueError("must be non-negative")

    @property
    def tidal_range(self):
        return self._tidal_range

    @tidal_range.setter
    def tidal_range(self, value):
        if value >= 0.0:
            self._tidal_range = float(value)
        else:
            raise ValueError("tidal range must be non-negative")

    @property
    def saturating_light(self):
        return self._saturating_light

    @saturating_light.setter
    def saturating_light(self, value):
        if value >= 0.0:
            self._saturating_light = float(value)
        else:
            raise ValueError("surface light must be non-negative")

    @property
    def sea_level(self):
        return self.grid.at_grid["sea_level__elevation"]

    @sea_level.setter
    def sea_level(self, sea_level):
        self.grid.at_grid["sea_level__elevation"] = float(sea_level)

    def _create_carbonate_thickness_field(self):
        """Create and return optional carbonate_thickness field."""
        return self.grid.add_zeros("carbonate_thickness", at="node")

    def calc_carbonate_production_rate(self):
        """Update carbonate production rate and store in field.

        Returns
        -------
        float array x number of grid nodes
            Reference to updated carbonate_production_rate field
        """
        self._depth[:] = self.sea_level - self._grid.at_node["topographic__elevation"]
        self._depth.clip(min=-2.0 * self._tidal_range - _EPSILON, out=self._depth)
        self._carb_prod_rate[self.grid.core_nodes] = (
            self._max_carbonate_production_rate
            * np.tanh(
                self.surface_light
                * np.exp(
                    -self.extinction_coefficient * self._depth[self.grid.core_nodes]
                )
                / self.saturating_light
            )
        )
        self._carb_prod_rate *= smooth_heaviside(self._depth, width=self.tidal_range)
        return self._carb_prod_rate

    def produce_carbonate(self, dt):
        """Grow carbonate for one time step & add to carbonate thickness field.

        If optional carbonate_thickness field does not already exist, this
        function creates it and initializes all values to zero.

        Returns
        -------
        float array x number of grid nodes
            Reference to updated carbonate_thickness field
        """
        # if self._carbonate_thickness is None:
        #     self._carbonate_thickness = self._create_carbonate_thickness_field()

        self.calc_carbonate_production_rate()
        # print(self._depth)
        # print(self._carb_prod_rate)
        added_thickness = self._carb_prod_rate * dt
        # print(added_thickness)
        self._carbonate_thickness += added_thickness
        # print(self._carbonate_thickness)
        self.grid.at_node["topographic__elevation"] += added_thickness
        # print(self.grid.at_node["topographic__elevation"])

        return self._carbonate_thickness

    def run_one_step(self, dt):
        """Advance by one time step.

        Simply calls the produce_carbonate() function.

        Parameters
        ----------
        dt : float
            Time step duration (normally in years)
        """
        self.produce_carbonate(dt)

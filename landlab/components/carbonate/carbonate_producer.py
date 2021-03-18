#! /usr/bin/env python
from contextlib import contextmanager

import numpy as np
from landlab import Component


@contextmanager
def set_numpy_err(*args, **kwds):
    settings = np.seterr(*args, **kwds)
    try:
        yield settings
    finally:
        np.seterr(**settings)


def smooth_heaviside(x, width=0.5, out=None):
    if width > 0.0:
        with set_numpy_err(over="ignore"):
            return np.divide(1.0, (1.0 + np.exp(-2.0 / width * np.asarray(x))), out=out)
    else:
        return np.heaviside(x, 0.5)


class CarbonateProducer(Component):
    """When I properly understand how to write the code, I'll write this intro ;-)"""

    _name = "The Producer of Carbonate"

    _time_units = "y"

    _info = {
        "sea_level__elevation": {
            "dtype": "float",
            "intent": "in",
            "optional": False,
            "units": "m",
            "mapping": "grid",
            "doc": "Position of sea level",
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
            "doc": "Carbonate production rate in latest time step",
        },
    }

    def __init__(
        self,
        grid,
        extinction_coefficient=1.0,
        max_carbonate_production_rate=1000.0,
        saturating_light=1.0,
        surface_light=1.0,
        tidal_range=0.0,
    ):
        """
        Parameters
        ----------
        grid: ModelGrid (RasterModelGrid, HexModelGrid, etc.)
            A landlab grid.
        """
        self.extinction_coefficient = extinction_coefficient
        self.max_carbonate_production_rate = max_carbonate_production_rate
        self.surface_light = surface_light
        self.tidal_range = tidal_range
        self.saturating_light = saturating_light

        super(CarbonateProducer, self).__init__(grid)

        if "carbonate_production_rate" not in self.grid.at_node:
            self.grid.add_empty("carbonate_production_rate", at="node")
        self.grid.at_node["carbonate_production_rate"].fill(0.0)

    @property
    def extinction_coefficient(self):
        return self._extinction_coefficient

    @extinction_coefficient.setter
    def extinction_coefficient(self, value):
        if value >= 0.0:
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
    def carbonate_thickness(self):
        return self._carb_thickness

    @property
    def sea_level(self):
        return self.grid.at_grid["sea_level__elevation"]

    @sea_level.setter
    def sea_level(self, sea_level):
        self.grid.at_grid["sea_level__elevation"] = sea_level

    def calc_production_with_depth(self, water_depth, out=None):
        """Return weighting factor for carbonate production.

        If there is no tidal range, then the weight factor is 1 at sea level,
        then decreases with increasing water depth following Bosscher and Schlager (1992),
        and 0 if above sea level. If there is a tidal range, then
        a tanh function is used to weight production across mean sea level, so
        that there is some degree of production for water depths within the
        tidal range (less above, more below). The nature of the tanh function
        is such that the production is about 95% of its maximum value at a depth
        of 1.5x the mean tidal range, and 5% of its maximum value at a height
        of 1.5x the mean tidal range above mean sea level.

        Parameters
        ----------
        water_depth : float array
            Depth of water relative to mean sea level (m) (can be negative)

            also needs, but not yet defined/coded
            surface_light
            extinction_coeff
            saturating_light

        Returns
        -------
        df : float array
            Weight factor ranging from 0 to 1.
        """
        water_depth = np.asarray(water_depth)
        if out is None:
            production = np.empty_like(water_depth)
        else:
            production = out

        under_water = water_depth >= 0.0
        production[under_water] = np.tanh(
            self.surface_light
            * np.exp(-self.extinction_coefficient * water_depth[under_water])
            / self.saturating_light
        )
        production[~under_water] = 0.0

        return production

    def calc_carbonate_production_at_node(self, out=None):
        """Calculate and store carbonate produced thickness, so strictly speaking accumulation.

        Returns
        -------
        carb_thick : float array
            produced carbonate thickness, m
        """
        water_depth = (
            self.sea_level
            + self.tidal_range
            - self._grid.at_node["topographic__elevation"]
        )
        fraction = self.calc_production_with_depth(
            water_depth, out=out
        ) * smooth_heaviside(water_depth, width=self.tidal_range)
        return np.multiply(self.max_carbonate_production_rate, fraction, out=out)

    def run_one_step(self):
        """Advance by one time step."""
        self.calc_carbonate_production_at_node(
            out=self._grid.at_node["carbonate_production_rate"]
        )

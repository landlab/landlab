#!/usr/bin/env python
"""Created on Fri Apr  8 08:32:48 2016.

@author: RCGlade
"""

import numpy as np

from landlab import Component


class ExponentialWeatherer(Component):
    """Calculate exponential weathering of bedrock on hillslopes.

    Uses exponential soil production function in the style of Ahnert (1976).

    Consider that ``w0`` is the maximum soil production rate and
    that ``w_star`` is the characteristic soil production depth. The
    soil production rate ``w0`` is given as a function of the soil
    depth::

        soil_production =  w0 * exp(-soil__depth / w_star)

    The :class:`~.ExponentialWeatherer` only calculates soil production at core nodes.

    An alternative version which uses the analytical integral of
    production through time is available at the component
    :py:class:`~.ExponentialWeathererIntegrated`.

    Examples
    --------
    >>> import numpy as np
    >>> from landlab import RasterModelGrid
    >>> from landlab.components import ExponentialWeatherer
    >>> mg = RasterModelGrid((5, 5))
    >>> soilz = mg.add_zeros("soil__depth", at="node")
    >>> soilrate = mg.add_ones("soil_production__rate", at="node")
    >>> expw = ExponentialWeatherer(mg)
    >>> expw.calc_soil_prod_rate()
    >>> np.allclose(mg.at_node["soil_production__rate"], 1.0)
    True

    References
    ----------
    **Required Software Citation(s) Specific to this Component**

    Barnhart, K., Glade, R., Shobe, C., Tucker, G. (2019). Terrainbento 1.0: a
    Python package for multi-model analysis in long-term drainage basin
    evolution. Geoscientific Model Development  12(4), 1267--1297.
    https://dx.doi.org/10.5194/gmd-12-1267-2019

    **Additional References**

    Ahnert, F. (1976). Brief description of a comprehensive three-dimensional
    process-response model of landform development Z. Geomorphol. Suppl.  25,
    29 - 49.

    Armstrong, A. (1976). A three dimensional simulation of slope forms.
    Zeitschrift für Geomorphologie  25, 20 - 28.
    """

    _name = "ExponentialWeatherer"

    _unit_agnostic = True

    _cite_as = """
    @article{barnhart2019terrain,
      author = {Barnhart, Katherine R and Glade, Rachel C and Shobe, Charles M
                and Tucker, Gregory E},
      title = {{Terrainbento 1.0: a Python package for multi-model analysis in
                long-term drainage basin evolution}},
      doi = {10.5194/gmd-12-1267-2019},
      pages = {1267---1297},
      number = {4},
      volume = {12},
      journal = {Geoscientific Model Development},
      year = {2019},
    }
    """

    _info = {
        "soil__depth": {
            "dtype": float,
            "intent": "in",
            "optional": False,
            "units": "m",
            "mapping": "node",
            "doc": "Depth of soil or weathered bedrock",
        },
        "soil_production__rate": {
            "dtype": float,
            "intent": "out",
            "optional": False,
            "units": "m/yr",
            "mapping": "node",
            "doc": "rate of soil production at nodes",
        },
    }

    def __init__(
        self, grid, soil_production_maximum_rate=1.0, soil_production_decay_depth=1.0
    ):
        """
        Parameters
        ----------
        grid: ModelGrid
            Landlab ModelGrid object
        soil_production_maximum_rate : float, array of float
            Maximum weathering rate for bare bedrock
        soil_production_decay_depth : float, array of float
            Characteristic weathering depth
        """
        super().__init__(grid)

        # Store grid and parameters
        # soil_production_maximum_rate, value set in setter below
        self._w0 = None
        # soil_production_decay_depth, value set in setter below
        self._wstar = None

        self.maximum_weathering_rate = soil_production_maximum_rate
        self.decay_depth = soil_production_decay_depth

        # weathering rate
        if "soil_production__rate" not in grid.at_node:
            grid.add_zeros("soil_production__rate", at="node")

    def calc_soil_prod_rate(self):
        """Calculate soil production rate."""
        core = self._grid.core_nodes
        self.grid.at_node["soil_production__rate"][core] = self._w0[core] * np.exp(
            -self.grid.at_node["soil__depth"][core] / self._wstar[core]
        )

    def run_one_step(self):
        """

        Parameters
        ----------
        dt: float
            Used only for compatibility with standard run_one_step.
        """
        self.calc_soil_prod_rate()

    @property
    def maximum_weathering_rate(self):
        """Maximum rate of weathering (m/yr)."""
        return self._w0

    @maximum_weathering_rate.setter
    def maximum_weathering_rate(self, new_val):
        if np.any(new_val <= 0):
            raise ValueError("Maximum weathering rate must be positive.")
        self._w0 = np.broadcast_to(new_val, self.grid.number_of_nodes)

    @property
    def decay_depth(self):
        """Maximum rate of weathering (m/yr)."""
        return self._wstar

    @decay_depth.setter
    def decay_depth(self, new_val):
        if np.any(new_val <= 0):
            raise ValueError("Maximum decay depth must be positive.")
        self._wstar = np.broadcast_to(new_val, self.grid.number_of_nodes)

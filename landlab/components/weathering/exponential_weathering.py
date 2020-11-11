#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""Created on Fri Apr  8 08:32:48 2016.

@author: RCGlade
"""

import numpy as np

from landlab import Component


class ExponentialWeatherer(Component):

    r"""
    This component implements exponential weathering of bedrock on hillslopes.
    Uses exponential soil production function in the style of Ahnert (1976).

    Consider that :math:`w_0` is the maximum soil production rate and
    that :math:`d^*` is the characteristic soil production depth. The
    soil production rate :math:`w` is given as a function of the soil
    depth :math:`d`,

    .. math::

        w = w_0 \exp{-\frac{d}{d^*}} \;.

    The `ExponentialWeatherer` only calculates soil production at core nodes.

    An alternative version which uses the analytical integral of
    production through time is available at the component
    :py:class:`ExponentialWeathererIntegrated <landlab.components.ExponentialWeathererIntegrated>`.

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
    >>> np.allclose(mg.at_node['soil_production__rate'], 1.)
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
      author = {Barnhart, Katherine R and Glade, Rachel C and Shobe, Charles M and Tucker, Gregory E},
      title = {{Terrainbento 1.0: a Python package for multi-model analysis in long-term drainage basin evolution}},
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
        self, grid, soil_production__maximum_rate=1.0, soil_production__decay_depth=1.0
    ):
        """
        Parameters
        ----------
        grid: ModelGrid
            Landlab ModelGrid object
        soil_production__maximum_rate : float
            Maximum weathering rate for bare bedrock
        soil_production__decay_depth : float
            Characteristic weathering depth
        """
        super().__init__(grid)

        # Store grid and parameters

        self._wstar = soil_production__decay_depth
        self._w0 = soil_production__maximum_rate

        # Create fields:
        # soil depth
        self._depth = grid.at_node["soil__depth"]

        # weathering rate
        if "soil_production__rate" in grid.at_node:
            self._soil_prod_rate = grid.at_node["soil_production__rate"]
        else:
            self._soil_prod_rate = grid.add_zeros("soil_production__rate", at="node")

    def calc_soil_prod_rate(self):
        """Calculate soil production rate."""
        # apply exponential function
        self._soil_prod_rate[self._grid.core_nodes] = self._w0 * np.exp(
            -self._depth[self._grid.core_nodes] / self._wstar
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
        if new_val <= 0:
            raise ValueError("Maximum weathering rate must be positive.")
        self._w0 = new_val

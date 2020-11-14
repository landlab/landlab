#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""Created on Fri Apr  8 08:32:48 2016.

@author: RCGlade
@author: dylanward
Integrated version created by D. Ward on Tue Oct 27 2020
"""

import numpy as np

from landlab import Component


class ExponentialWeathererIntegrated(Component):

    r"""
    This component implements exponential weathering of bedrock on
    hillslopes. Uses exponential soil production function in the style
    of Ahnert (1976).

    Consider that :math:`w_0` is the maximum soil production rate and
    that :math:`d^*` is the characteristic soil production depth. The
    soil production rate :math:`w` is given as a function of the soil
    depth :math:`d`,

    .. math::

        w = w_0 \exp{-\frac{d}{d^*}} \;.

    The `ExponentialWeathererIntegrated` uses the analytical solution
    for the amount of soil produced by an exponential weathering
    function over a timestep dt, and returns both the thickness of
    bedrock weathered and the thickness of soil produced. The solution
    accounts for the reduction in rate over the timestep due to the
    increasing depth. This enables accuracy over arbitrarily large
    timesteps, and better compatiblity with the `run_one_step()`
    interface.

    Compared to 'ExponentialWeatherer', upon which it is based...

    - This maintains the field I/O behavior of the original, but adds
      new return fields for the weathered thickness and soil produced
      thickness.
    - Density adjustments are needed inside the integral and the
      density ratio is intialized on instantiation. The default value
      of 1.0 assumes no change in density.
    - Returns both weathered depth of bedrock and produced depth of
      soil over the timestep.
    - The primary `soil__depth` field that is input is NOT updated by
      the component.

    This is left as an exercise for the model driver, as different
    applications may want to integrate soil depth and weathering in
    different sequences among other processes.

    - SHOULD maintain drop-in compatiblity with the plain
      :py:class:`ExponentialWeatherer <landlab.components.ExponentialWeatherer>`,
      just import and instantiate this one instead and existing code
      should work with no side effects other than the creation of the
      two additional (zeros) output fields.

    Examples
    --------
    >>> import numpy as np
    >>> from landlab import RasterModelGrid
    >>> from landlab.components import ExponentialWeathererIntegrated
    >>> mg = RasterModelGrid((5, 5))
    >>> soilz = mg.add_zeros("soil__depth", at="node")
    >>> soilrate = mg.add_ones("soil_production__rate", at="node")
    >>> expw = ExponentialWeathererIntegrated(mg)
    >>> dt = 1000
    >>> expw.run_one_step(dt)
    >>> np.allclose(mg.at_node['soil_production__rate'][mg.core_nodes], 1.)
    True
    >>> np.allclose(mg.at_node['soil_production__dt_produced_depth'][mg.core_nodes], 6.9088)
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

    _name = "ExponentialWeathererIntegrated"

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
        "soil_production__dt_produced_depth": {
            "dtype": float,
            "intent": "out",
            "optional": False,
            "units": "m",
            "mapping": "node",
            "doc": "thickness of soil produced at nodes over time dt",
        },
        "soil_production__dt_weathered_depth": {
            "dtype": float,
            "intent": "out",
            "optional": False,
            "units": "m",
            "mapping": "node",
            "doc": "thickness of bedrock weathered at nodes over time dt",
        },
    }

    def __init__(
        self,
        grid,
        soil_production__maximum_rate=1.0,
        soil_production__decay_depth=1.0,
        soil_production__expansion_factor=1.0,
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
        soil_production__expansion_factor : float
            Expansion ratio of regolith (from relative densities of
            rock and soil)
        """
        super().__init__(grid)

        # Store grid and parameters

        self._wstar = soil_production__decay_depth
        self._w0 = soil_production__maximum_rate
        self._fexp = soil_production__expansion_factor

        # Create fields:
        # soil depth
        self._depth = grid.at_node["soil__depth"]

        # weathering rate
        if "soil_production__rate" in grid.at_node:
            self._soil_prod_rate = grid.at_node["soil_production__rate"]
        else:
            self._soil_prod_rate = grid.add_zeros("soil_production__rate", at="node")

        # soil produced total over dt
        if "soil_production__dt_produced_depth" in grid.at_node:
            self._soil_prod_total = grid.at_node["soil_production__dt_produced_depth"]
        else:
            self._soil_prod_total = grid.add_zeros(
                "soil_production__dt_produced_depth", at="node"
            )

        # bedrock weathering total over dt
        if "soil_production__dt_weathered_depth" in grid.at_node:
            self._rock_weathered_total = grid.at_node[
                "soil_production__dt_weathered_depth"
            ]
        else:
            self._rock_weathered_total = grid.add_zeros(
                "soil_production__dt_weathered_depth", at="node"
            )

    def calc_soil_prod_rate(self):
        """Calculate soil production rate."""
        # apply exponential function
        self._soil_prod_rate[self._grid.core_nodes] = self._w0 * np.exp(
            -self._depth[self._grid.core_nodes] / self._wstar
        )

    def _calc_dt_production_total(self, dt):
        """Calculate integrated production over 1 timestep dt"""
        # analytical solution
        self._soil_prod_total[self._grid.core_nodes] = self._wstar * np.log(
            (
                self._fexp
                * self._soil_prod_rate[self._grid.core_nodes]
                * dt
                / self._wstar
            )
            + 1
        )
        # and back-convert to find rock thickness converted over the timestep:
        self._rock_weathered_total[self._grid.core_nodes] = (
            self._soil_prod_total[self._grid.core_nodes] / self._fexp
        )

    def run_one_step(self, dt=0):
        """
        Parameters
        ----------
        dt: float
            Used only for compatibility with standard run_one_step.
            If dt is not provided, the default of zero maintains backward compatibility
        """
        self.calc_soil_prod_rate()
        self._calc_dt_production_total(dt)

    @property
    def maximum_weathering_rate(self):
        """Maximum rate of weathering (m/yr)."""
        return self._w0

    @maximum_weathering_rate.setter
    def maximum_weathering_rate(self, new_val):
        if new_val <= 0:
            raise ValueError("Maximum weathering rate must be positive.")
        self._w0 = new_val

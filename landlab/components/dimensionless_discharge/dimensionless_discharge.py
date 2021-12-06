#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Calculate dimensionless discharge of stream sections based on Tang et
al. (2019)
"""

import numpy as np

from landlab import Component


class DimensionlessDischarge(Component):
    r"""Component that calculates dimensionless discharge of stream
     segments.

    The dimensionless discharge model calculates the unitless  discharge
    value of streams and can be used to help determine locations of
    debris flows. It uses an equation from Tang et al. (2019) to
    calculate the dimensionless discharge as well as the threshold for
    whether a debris flow will occur for a specified location.

    .. math::

        q* = \frac{q}{\sqrt{\frac{\rho_s-\rho}{\rho}*g*D_{50}^3}}

     where :math:`q*` is the dimensionless discharge value for a stream
     segment, :math:`q` is flux in a stream segment,  :math:`\rho_s`
     is soil density, :math:`\rho` is water density, :math:`g` is the
     gravitational constant, and :math:`D_50` is the average sediment
     partical size in the stream segment area.

     Constants C and N are coefficients used in the slope-dependent
     equation

     .. math::

        q*_{thresold} = \frac{C}{(tan(\theta))^N}

     to determine whether the dimensionless discharge calculated
     exceeds thresholds for a sediment-laden (Upper Limit) or
     water-producing only (Lower Limit) debris flow.  C and N are
     empirically-derived by comparing dimensionless discharge estimates
     against previous debris flow events.  The values will vary based on
     local geology and soil type. In southern California values of C=12
     and N=0.85 (upper limits, Tang et al, 2019) and C=4.29, N=0.78
     (lower limits, Tang et al, 2019) have been used, while values of
     C=0.195 and N=1.27 have been in used in the Italian Dolomites
     (Gregoretti and Fontana, 2008). Default values are C=12 and N=0.85.

     Parameters
     ----------
     soil_density : float list
         Density of soil in watershed (kg/m^3)
     water_density : float
         density of water in waterhed (kg/m^3)
     C : float
         Numerator of the debris flow threshold equation; Empirically
         derived constant (See Tang et al. 2019)
     N : float
         Exponent for slope in the denominator of the debris flow
         threshold equation; Empirically derived constant (See Tang et al. 2019)

     Examples
     --------
     >>> from landlab.components import DimensionlessDischarge
     >>> from landlab import RasterModelGrid
     >>> import random
     >>> watershed_grid = RasterModelGrid((3, 3))
     >>> flux = watershed_grid.add_ones('node', 'flux')
     >>> d50 = watershed_grid.add_ones('node', 'd50')
     >>> watershed_grid.at_node['dem_values'] = np.array([[1.1, 2, 3, 4, 2, 3, 4, 5, 3]])
     >>> dd = DimensionlessDischarge(watershed_grid)
     >>> dd.run_one_step()
     >>> print(watershed_grid.at_node['dimensionless_discharge'])
     [ 0.55372743  0.55372743  0.55372743  0.55372743  0.55372743
             0.55372743  0.55372743  0.55372743  0.55372743]



     References
     ----------
     Tang, H., McGuire, L. A., Rengers, F. K., Kean, J. W., Staley,
     D. M., & Smith, J. B. (2019). Developing and Testing Physically
     Based Triggering Thresholds for Runoff-Generated Debris Flows.
     Geophysical Research Letters, 46(15), 8830–8839.
     https://doi.org/10.1029/2019GL083623

    """

    _name = "DimensionlessDischargeModel"

    _unit_agnostic = False

    _info = {
        "dimensionless_discharge": {
            "dtype": float,
            "intent": "out",
            "optional": False,
            "units": "none",
            "mapping": "node",
            "doc": "Dimensionless discharge value for a stream segment.",
        },
        "dimensionless_discharge_above_threshold": {
            "dtype": bool,
            "intent": "out",
            "optional": False,
            "units": "none",
            "mapping": "node",
            "doc": "Dimensionless discharge value for a stream segment.",
        },
        "flux": {
            "dtype": float,
            "intent": "in",
            "optional": False,
            "units": "m/s",
            "mapping": "node",
            "doc": "flux in each stream segment",
        },
        "d50": {
            "dtype": float,
            "intent": "in",
            "optional": False,
            "units": "m",
            "mapping": "node",
            "doc": "soil grain size average in stream segment",
        },
        "dem_values": {
            "dtype": float,
            "intent": "in",
            "optional": False,
            "units": "m",
            "mapping": "node",
            "doc": "elevation stream segment",
        },
    }

    def __init__(self, grid, soil_density=1330, water_density=997.9, C=12.0, N=0.85):
        """Initialize the DimensionlessDischarge.

        Parameters
        ----------
        soil_density : float, optional (defaults to 1330)
            density of soil (kg/m^3)
        water_density : float, optional (defaults to 997.9)
            density of water (kg/m^3)
        C : float, optional (defaults to 12.0)
            Numerator of the debris flow threshold equation; Empirically
            derived constant (See Tang et al. 2019)
        N : float, optional (defaults to 0.85)
            Exponent for slope in the denominator of the debris flow
            threshold equation; Empirically derived constant (See Tang
            et al. 2019)
        """

        super().__init__(grid)

        # Store parameters
        self._soil_density = soil_density
        self._C = [C] * grid.number_of_nodes
        self._N = [N] * grid.number_of_nodes
        # change DEM values into slope of a stream segment
        self._stream_slopes = grid.calc_slope_at_node(elevs="dem_values")
        self.water_density = water_density
        self.gravity = 9.8

        # set threshold values for each segment
        _ = self.grid.add_zeros("node", "dimensionless_discharge")
        _ = self.grid.add_zeros("node", "dimensionless_discharge_above_threshold")
        self.grid.at_node["dimensionless_discharge_above_threshold"] = np.array(
            [[False] * self.grid.number_of_nodes]
        )
        _ = self.grid.add_zeros("node", "dimensionless_discharge_threshold_value")
        self.grid.at_node["dimensionless_discharge_threshold_value"] = self._C / (
            self._stream_slopes ** self._N
        )

    def run_one_step(self):

        self.grid.at_node["dimensionless_discharge"] = self.grid.at_node[
            "flux"
        ] / np.sqrt(
            ((self._soil_density - self.water_density) / self.water_density)
            * self.gravity
            * (self.grid.at_node["d50"] ** 3)
        )

        self.grid.at_node["dimensionless_discharge_above_threshold"] = [
            True
            if self.grid.at_node["dimensionless_discharge"][i]
            >= self.grid.at_node["dimensionless_discharge_threshold_value"][i]
            else False
            for i in range(self.grid.number_of_nodes)
        ]

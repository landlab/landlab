#!/usr/bin/env python3
"""
Calculate dimensionless discharge of stream sections based on Tang et
al. (2019)
"""

import numpy as np
from scipy import constants

from landlab import Component
from landlab.utils.return_array import return_array_at_node


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

     Examples
     --------
     >>> from landlab.components import DimensionlessDischarge
     >>> from landlab import RasterModelGrid
     >>> import random
     >>> watershed_grid = RasterModelGrid((3, 3))
     >>> surface_water__unit_discharge = watershed_grid.add_ones(
     ...     "surface_water__unit_discharge", at="node"
     ... )
     >>> d50 = watershed_grid.add_ones(
     ...     "channel_bottom_sediment_grain__d50_diameter", at="node"
     ... )
     >>> watershed_grid.at_node["topographic__elevation"] = np.array(
     ...     [[1.1, 2, 3, 4, 2, 3, 4, 5, 3]]
     ... )
     >>> dd = DimensionlessDischarge(watershed_grid, gravity=9.8)
     >>> dd.run_one_step()
     >>> watershed_grid.at_node["dimensionless_discharge"]
     array([0.55372743, 0.55372743, 0.55372743, 0.55372743, 0.55372743,
            0.55372743, 0.55372743, 0.55372743, 0.55372743])

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
            "doc": (
                "True if dimensionless discharge value is above threshold value, "
                "false otherwise."
            ),
        },
        "dimensionless_discharge_threshold": {
            "dtype": float,
            "intent": "out",
            "optional": False,
            "units": "none",
            "mapping": "node",
            "doc": "Dimensionless discharge threshold for each stream segment.",
        },
        "surface_water__unit_discharge": {
            "dtype": float,
            "intent": "in",
            "optional": False,
            "units": "m**2/s",
            "mapping": "node",
            "doc": "Volumetric discharge of surface water per unit width",
        },
        "channel_bottom_sediment_grain__d50_diameter": {
            "dtype": float,
            "intent": "in",
            "optional": True,
            "units": "m",
            "mapping": "node",
            "doc": "soil grain size average in stream segment",
        },
        "topographic__elevation": {
            "dtype": float,
            "intent": "inout",
            "optional": False,
            "units": "m",
            "mapping": "node",
            "doc": "Land surface topographic elevation",
        },
    }

    def __init__(
        self,
        grid,
        soil_density=1330,
        water_density=997.9,
        C=12.0,
        N=0.85,
        gravity=constants.g,
        channel_bottom_sediment_grain__d50_diameter="channel_bottom_sediment_grain__d50_diameter",
    ):
        """Initialize the DimensionlessDischarge.

        Parameters
        ----------
        soil_density : float, field name, or array, optional
            density of soil (kg/m^3) (defaults to 1330)
        water_density : float, optional
            density of water (kg/m^3) (defaults to 997.9)
        C : float, optional
            Numerator of the debris flow threshold equation; Empirically
            derived constant (See Tang et al. 2019). (defaults to 12.0)
        N : float, optional
            Exponent for slope in the denominator of the debris flow
            threshold equation; Empirically derived constant (See Tang
            et al. 2019). (defaults to 0.85)
        gravity : float, optional
            (defaults to scipy's standard acceleration of gravity)
        channel_bottom_sediment_grain__d50_diameter : float, field_name, or array
            defaults to the field name
            "channel_bottom_sediment_grain__d50_diameter"
        """

        super().__init__(grid)

        # Store parameters

        # scalar parameters.
        self._c = C
        self._n = N
        self._water_density = water_density
        self._gravity = gravity

        # get topographic__elevation values and change into slope of a
        # stream segment
        self._stream_slopes = self._elevationToSlope()

        # both soil density and d50 can be float, array, or field name.
        self._soil_density = return_array_at_node(self.grid, soil_density)
        self._channel_bottom_sediment_grain__d50_diameter = return_array_at_node(
            self.grid, channel_bottom_sediment_grain__d50_diameter
        )

        # create output fields
        self.grid.add_zeros("dimensionless_discharge", at="node")
        self.grid.add_full(
            "dimensionless_discharge_above_threshold", False, at="node", dtype=bool
        )
        self.grid.add_zeros("dimensionless_discharge_threshold", at="node", dtype=float)

        # calculate threshold values for each segment
        self._calc_threshold()

    def _calc_threshold(self):
        self.grid.at_node["dimensionless_discharge_threshold"] = self._c / (
            self._stream_slopes**self._n
        )

    def _elevationToSlope(self):
        return self.grid.calc_slope_at_node(elevs="topographic__elevation")

    def run_one_step(self):
        # update slopes
        self._stream_slopes = self._elevationToSlope()

        # recalculate threshold with new slopes.
        self._calc_threshold()

        self.grid.at_node["dimensionless_discharge"] = self.grid.at_node[
            "surface_water__unit_discharge"
        ] / np.sqrt(
            ((self._soil_density - self._water_density) / self._water_density)
            * self._gravity
            * (self._channel_bottom_sediment_grain__d50_diameter**3)
        )

        self.grid.at_node["dimensionless_discharge_above_threshold"] = (
            self.grid.at_node["dimensionless_discharge"]
            > self.grid.at_node["dimensionless_discharge_threshold"]
        )

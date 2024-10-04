#!/usr/bin/env python3
"""
Model bedrock incision and gravel transport and abrasion in a network of rivers.

@author: gtucker
"""

import numpy as np

from landlab import Component
from landlab import HexModelGrid
from landlab.grid.diagonals import DiagonalsMixIn

_DT_MAX = 1.0e-2
_ONE_SIXTH = 1.0 / 6.0
_SEVEN_SIXTHS = 7.0 / 6.0


class GravelBedrockEroder(Component):
    """Drainage network evolution of rivers with gravel alluvium overlying bedrock.

    Model drainage network evolution for a network of rivers that have
    a layer of gravel alluvium overlying bedrock.

    :class:`~.GravelBedrockEroder` is designed to operate together with a flow-routing
    component such as :class:`~.FlowAccumulator`, so that each grid node has
    a defined flow direction toward one of its neighbor nodes. Each core node
    is assumed to contain one outgoing fluvial channel, and (depending on
    the drainage structure) zero, one, or more incoming channels. These channels are
    treated as effectively sub-grid-scale features that are embedded in valleys
    that have a width of one grid cell.

    As with the :class:`~.GravelRiverTransporter` component, the rate of gravel
    transport out of a given node is calculated as the product of bankfull discharge,
    channel gradient (to the 7/6 power), a dimensionless transport coefficient, and
    an intermittency factor that represents the fraction of time that bankfull
    flow occurs. The derivation of the transport law is given by Wickert &
    Schildgen (2019), and it derives from the assumption that channels are
    gravel-bedded and that they "instantaneously" adjust their width such that
    bankfull bed shear stress is just slightly higher than the threshold for
    grain motion. The substrate is assumed to consist entirely of gravel-size
    material with a given bulk porosity. The component calculates the loss of
    gravel-sized material to abrasion (i.e., conversion to finer sediment, which
    is not explicitly tracked) as a function of the volumetric transport rate,
    an abrasion coefficient with units of inverse length, and the local transport
    distance (for example, if a grid node is carrying a gravel load ``Qs`` to a
    neighboring node ``dx`` meters downstream, the rate of gravel loss in volume per
    time per area at the node will be ``beta * Qs * dx``, where ``beta`` is the abrasion
    coefficient).

    Sediment mass conservation is calculated across each entire
    grid cell. For example, if a cell has surface area ``A``, a total volume influx
    ``Qin``, and downstream transport rate ``Qs``, the resulting rate of change of
    alluvium thickness will be ``(Qin - Qs / (A * (1 - phi))``, plus gravel produced by
    plucking erosion of bedrock (``phi`` is porosity).

    Bedrock is eroded by a combination of abrasion and plucking. Abrasion per unit
    channel length is calculated as the product of volumetric sediment discharge
    and an abrasion coefficient. Sediment produced by abrasion is assumed to
    go into wash load that is removed from the model domain. Plucking is calculated
    using a discharge-slope expression, and a user-defined fraction of plucked
    material is added to the coarse alluvium.

    Parameters
    ----------
    grid : ModelGrid
        A Landlab model grid object
    intermittency_factor : float (default 0.01)
        Fraction of time that bankfull flow occurs
    transport_coefficient : float (default 0.041)
        Dimensionless transport efficiency factor (see Wickert & Schildgen 2019)
    abrasion_coefficient : float (default 0.0 1/m)
        Abrasion coefficient with units of inverse length
    sediment_porosity : float (default 0.35)
        Bulk porosity of bed sediment
    depth_decay_scale : float (default 1.0)
        Scale for depth decay in bedrock exposure function
    plucking_coefficient : float or (n_core_nodes,) array of float (default 1.0e-4 1/m)
        Rate coefficient for bedrock erosion by plucking
    coarse_fraction_from_plucking : float or (n_core_nodes,) array of float (default 1.0)
        Fraction of plucked material that becomes part of gravel sediment load

    Examples
    --------
    >>> from landlab import RasterModelGrid
    >>> from landlab.components import FlowAccumulator
    >>> grid = RasterModelGrid((3, 3), xy_spacing=1000.0)
    >>> elev = grid.add_zeros("topographic__elevation", at="node")
    >>> sed = grid.add_zeros("soil__depth", at="node")
    >>> sed[4] = 300.0
    >>> grid.status_at_node[grid.perimeter_nodes] = grid.BC_NODE_IS_CLOSED
    >>> grid.status_at_node[5] = grid.BC_NODE_IS_FIXED_VALUE
    >>> fa = FlowAccumulator(grid, runoff_rate=10.0)
    >>> fa.run_one_step()
    >>> eroder = GravelBedrockEroder(grid, abrasion_coefficient=0.0005)
    >>> rock_elev = grid.at_node["bedrock__elevation"]
    >>> for _ in range(200):
    ...     rock_elev[grid.core_nodes] += 1.0
    ...     elev[grid.core_nodes] += 1.0
    ...     fa.run_one_step()
    ...     eroder.run_one_step(10000.0)
    ...
    >>> int(elev[4] * 100)
    2266
    """

    _name = "GravelBedrockEroder"

    _unit_agnostic = True

    _info = {
        "bedload_sediment__rate_of_loss_to_abrasion": {
            "dtype": float,
            "intent": "out",
            "optional": False,
            "units": "m/y",
            "mapping": "node",
            "doc": "Rate of bedload sediment volume loss to abrasion per unit area",
        },
        "bedload_sediment__volume_influx": {
            "dtype": float,
            "intent": "out",
            "optional": False,
            "units": "m**3/y",
            "mapping": "node",
            "doc": "Volumetric incoming streamwise bedload sediment transport rate",
        },
        "bedload_sediment__volume_outflux": {
            "dtype": float,
            "intent": "out",
            "optional": False,
            "units": "m**3/y",
            "mapping": "node",
            "doc": "Volumetric outgoing streamwise bedload sediment transport rate",
        },
        "bedrock__abrasion_rate": {
            "dtype": float,
            "intent": "out",
            "optional": False,
            "units": "m/y",
            "mapping": "node",
            "doc": "rate of bedrock lowering by abrasion",
        },
        "bedrock__elevation": {
            "dtype": float,
            "intent": "out",
            "optional": False,
            "units": "m",
            "mapping": "node",
            "doc": "elevation of the bedrock surface",
        },
        "bedrock__exposure_fraction": {
            "dtype": float,
            "intent": "out",
            "optional": False,
            "units": "-",
            "mapping": "node",
            "doc": "fractional exposure of bedrock",
        },
        "bedrock__plucking_rate": {
            "dtype": float,
            "intent": "out",
            "optional": False,
            "units": "m/y",
            "mapping": "node",
            "doc": "rate of bedrock lowering by plucking",
        },
        "bedrock__lowering_rate": {
            "dtype": float,
            "intent": "out",
            "optional": False,
            "units": "m/y",
            "mapping": "node",
            "doc": "Rate of lowering of bedrock surface",
        },
        "flow__link_to_receiver_node": {
            "dtype": int,
            "intent": "in",
            "optional": False,
            "units": "-",
            "mapping": "node",
            "doc": "ID of link downstream of each node, which carries the discharge",
        },
        "flow__receiver_node": {
            "dtype": int,
            "intent": "in",
            "optional": False,
            "units": "-",
            "mapping": "node",
            "doc": "Node array of receivers (node that receives flow from current node)",
        },
        "flow__upstream_node_order": {
            "dtype": int,
            "intent": "in",
            "optional": False,
            "units": "-",
            "mapping": "node",
            "doc": "Node array containing downstream-to-upstream ordered list of node IDs",
        },
        "sediment__rate_of_change": {
            "dtype": float,
            "intent": "out",
            "optional": False,
            "units": "m/y",
            "mapping": "node",
            "doc": "Time rate of change of sediment thickness",
        },
        "soil__depth": {
            "dtype": float,
            "intent": "in",
            "optional": False,
            "units": "m",
            "mapping": "node",
            "doc": "Depth of soil or weathered bedrock",
        },
        "surface_water__discharge": {
            "dtype": float,
            "intent": "in",
            "optional": False,
            "units": "m**3/y",
            "mapping": "node",
            "doc": "Volumetric discharge of surface water",
        },
        "topographic__elevation": {
            "dtype": float,
            "intent": "inout",
            "optional": False,
            "units": "m",
            "mapping": "node",
            "doc": "Land surface topographic elevation",
        },
        "topographic__steepest_slope": {
            "dtype": float,
            "intent": "in",
            "optional": False,
            "units": "-",
            "mapping": "node",
            "doc": "The steepest *downhill* slope",
        },
    }

    def __init__(
        self,
        grid,
        intermittency_factor=0.01,
        transport_coefficient=0.041,
        abrasion_coefficient=0.0,
        sediment_porosity=0.35,
        depth_decay_scale=1.0,
        plucking_coefficient=1.0e-4,
        coarse_fraction_from_plucking=1.0,
    ):
        """Initialize GravelBedrockEroder."""

        super().__init__(grid)

        # Parameters
        self._trans_coef = transport_coefficient
        self._intermittency_factor = intermittency_factor
        self._abrasion_coef = abrasion_coefficient
        self._porosity_factor = 1.0 / (1.0 - sediment_porosity)
        self._depth_decay_scale = depth_decay_scale
        if (
            isinstance(plucking_coefficient, np.ndarray)
            and len(plucking_coefficient) == self.grid.number_of_nodes
        ):
            plucking_coefficient = plucking_coefficient[self.grid.core_nodes]
        self._plucking_coef = plucking_coefficient
        if (
            isinstance(coarse_fraction_from_plucking, np.ndarray)
            and len(coarse_fraction_from_plucking) == self.grid.number_of_nodes
        ):
            coarse_fraction_from_plucking = coarse_fraction_from_plucking[
                self.grid.core_nodes
            ]
        self._pluck_coarse_frac = coarse_fraction_from_plucking

        # Fields and arrays
        self._elev = grid.at_node["topographic__elevation"]
        self._sed = grid.at_node["soil__depth"]
        if "bedrock__elevation" in grid.at_node:
            self._bedrock__elevation = grid.at_node["bedrock__elevation"]
        else:
            self._bedrock__elevation = grid.add_zeros(
                "bedrock__elevation", at="node", dtype=float
            )
            self._bedrock__elevation[:] = self._elev - self._sed
        self._discharge = grid.at_node["surface_water__discharge"]
        self._slope = grid.at_node["topographic__steepest_slope"]
        self._receiver_node = grid.at_node["flow__receiver_node"]
        self._receiver_link = grid.at_node["flow__link_to_receiver_node"]
        super().initialize_output_fields()
        self._sediment_influx = grid.at_node["bedload_sediment__volume_influx"]
        self._sediment_outflux = grid.at_node["bedload_sediment__volume_outflux"]
        self._dHdt = grid.at_node["sediment__rate_of_change"]
        self._rock_lowering_rate = grid.at_node["bedrock__lowering_rate"]
        self._abrasion = grid.at_node["bedload_sediment__rate_of_loss_to_abrasion"]
        self._rock_exposure_fraction = grid.at_node["bedrock__exposure_fraction"]
        self._rock_abrasion_rate = grid.at_node["bedrock__abrasion_rate"]
        self._pluck_rate = grid.at_node["bedrock__plucking_rate"]

        self._setup_length_of_flow_link()

    def _setup_length_of_flow_link(self):
        """Set up a float or array containing length of the flow link from
        each node, which is needed for the abrasion rate calculations.
        """
        if isinstance(self.grid, HexModelGrid):
            self._flow_link_length_over_cell_area = (
                self.grid.spacing / self.grid.area_of_cell[0]
            )
            self._flow_length_is_variable = False
        elif isinstance(self.grid, DiagonalsMixIn):
            self._flow_length_is_variable = True
            self._grid_has_diagonals = True
            self._update_flow_link_length_over_cell_area()
        else:
            self._flow_length_is_variable = True
            self._grid_has_diagonals = False
            self._update_flow_link_length_over_cell_area()

    def _update_flow_link_length_over_cell_area(self):
        """Update the ratio of the length of link along which water flows out of
        each node to the area of the node's cell."""
        if self._grid_has_diagonals:
            flow_link_len = self.grid.length_of_d8
        else:
            flow_link_len = self.grid.length_of_link
        self._flow_link_length_over_cell_area = (
            flow_link_len[self._receiver_link[self.grid.core_nodes]]
            / self.grid.area_of_cell[self.grid.cell_at_node[self.grid.core_nodes]]
        )

    def calc_implied_depth(self, grain_diameter=0.01):
        """Utility function that calculates and returns water depth implied by
        slope and grain diameter, using Wickert & Schildgen (2019) equation 8.

        The equation is::

            h = ((rho_s - rho / rho)) * (1 + epsilon) * tau_c * (D / S)

        where the factors on the right are sediment and water density, excess
        shear-stress factor, critical Shields stress, grain diameter, and slope
        gradient. Here the prefactor on ``D/S`` assumes sediment density of 2650 kg/m3,
        water density of 1000 kg/m3, shear-stress factor of 0.2, and critical
        Shields stress of 0.0495, giving a value of 0.09801.

        Examples
        --------
        >>> from landlab import RasterModelGrid
        >>> from landlab.components import FlowAccumulator
        >>> grid = RasterModelGrid((3, 3), xy_spacing=1000.0)
        >>> elev = grid.add_zeros("topographic__elevation", at="node")
        >>> elev[3:] = 10.0
        >>> sed = grid.add_zeros("soil__depth", at="node")
        >>> sed[3:] = 100.0
        >>> fa = FlowAccumulator(grid)
        >>> fa.run_one_step()
        >>> eroder = GravelBedrockEroder(grid)
        >>> water_depth = eroder.calc_implied_depth()
        >>> int(water_depth[4] * 1000)
        98
        """
        depth_factor = 0.09801
        depth = np.zeros(self._grid.number_of_nodes)
        nonzero_slope = self._slope > 0.0
        depth[nonzero_slope] = (
            depth_factor * grain_diameter / self._slope[nonzero_slope]
        )
        return depth

    def calc_implied_width(self, grain_diameter=0.01, time_unit="y"):
        """Utility function that calculates and returns channel width implied by
        discharge, slope, and grain diameter, using Wickert & Schildgen (2019)
        equation 16.

        The equation is::

            b = kb * Q * S**(7/6) / D**(3/2)

        where the dimensional prefactor, which includes sediment and water
        density, gravitational acceleration, critical Shields stress, and the
        transport factor epsilon, is::

            kb = 0.17 g**(-1/2) (((rho_s - rho) / rho) (1 + eps) tau_c*)**(-5/3)

        Using ``g = 9.8 m/s2``, ``rho_s = 2650`` (quartz), ``rho = 1000 kg/m3``, ``eps = 0.2``,
        and ``tau_c* = 0.0495``, ``kb ~ 2.61 s/m**(1/2)``. Converting to years,
        ``kb = 8.26e-8``.

        Examples
        --------
        >>> from landlab import RasterModelGrid
        >>> from landlab.components import FlowAccumulator
        >>> grid = RasterModelGrid((3, 3), xy_spacing=10000.0)
        >>> elev = grid.add_zeros("topographic__elevation", at="node")
        >>> elev[3:] = 100.0
        >>> sed = grid.add_zeros("soil__depth", at="node")
        >>> sed[3:] = 100.0
        >>> fa = FlowAccumulator(grid)
        >>> fa.run_one_step()
        >>> eroder = GravelBedrockEroder(grid)
        >>> chan_width = eroder.calc_implied_width()
        >>> int(chan_width[4] * 100)
        3833
        >>> grid.at_node["surface_water__discharge"] *= 1.0 / (3600 * 24 * 365.25)
        >>> chan_width = eroder.calc_implied_width(time_unit="s")
        >>> int(chan_width[4] * 100)
        3838
        """
        if time_unit[0] == "y":
            width_fac = 8.26e-8
        else:
            width_fac = 2.61  # assume seconds if not years
        width = (
            width_fac
            * self._discharge
            * self._slope ** (7.0 / 6.0)
            / (grain_diameter**1.5)
        )
        return width

    def calc_rock_exposure_fraction(self):
        """Update the bedrock exposure fraction.

        The result is stored in the ``bedrock__exposure_fraction`` field.

        Examples
        --------
        >>> from landlab import RasterModelGrid
        >>> from landlab.components import FlowAccumulator
        >>> grid = RasterModelGrid((3, 4), xy_spacing=100.0)
        >>> elev = grid.add_zeros("topographic__elevation", at="node")
        >>> sed = grid.add_zeros("soil__depth", at="node")
        >>> sed[4] = 1000.0
        >>> sed[5] = 0.0
        >>> fa = FlowAccumulator(grid)
        >>> fa.run_one_step()
        >>> eroder = GravelBedrockEroder(grid)
        >>> eroder.calc_rock_exposure_fraction()
        >>> eroder._rock_exposure_fraction[4:6]
        array([0., 1.])
        >>> sed[4] = 1.0  # exposure frac should be 1/e ~ 0.3679
        >>> sed[5] = 2.0  # exposure frac should be 1/e^2 ~ 0.1353
        >>> eroder.calc_rock_exposure_fraction()
        >>> np.round(eroder._rock_exposure_fraction[4:6], 4)
        array([0.3679, 0.1353])
        """
        self._rock_exposure_fraction[:] = np.exp(-self._sed / self._depth_decay_scale)

    def calc_transport_rate(self):
        """Calculate and return bed-load transport rate.

        Calculation uses Wickert-Schildgen approach, and provides
        volume per time rate. Transport rate is modulated by available
        sediment, using the exponential function ``(1 - exp(-H / Hs))``,
        so that transport rate approaches zero as sediment thickness
        approaches zero. Rate is a volume per time. The result is
        stored in the ``bedload_sediment__volume_outflux`` field.

        Examples
        --------
        >>> from landlab import RasterModelGrid
        >>> from landlab.components import FlowAccumulator
        >>> grid = RasterModelGrid((3, 3), xy_spacing=100.0)
        >>> elev = grid.add_zeros("topographic__elevation", at="node")
        >>> elev[3:] = 1.0
        >>> sed = grid.add_zeros("soil__depth", at="node")
        >>> sed[3:] = 100.0
        >>> fa = FlowAccumulator(grid)
        >>> fa.run_one_step()
        >>> eroder = GravelBedrockEroder(grid)
        >>> eroder.calc_transport_rate()
        >>> round(eroder._sediment_outflux[4], 4)
        0.019
        """
        self._sediment_outflux[:] = (
            self._trans_coef
            * self._intermittency_factor
            * self._discharge
            * self._slope**_SEVEN_SIXTHS
            * (1.0 - self._rock_exposure_fraction)
        )

    def calc_abrasion_rate(self):
        """Update the volume rate of bedload loss to abrasion, per unit area.

        Here we use the average of incoming and outgoing sediment flux to
        calculate the loss rate to abrasion. The result is stored in the
        ``bedload_sediment__rate_of_loss_to_abrasion`` field.

        The factor dx (node spacing) appears in the denominator to represent
        flow segment length (i.e., length of the link along which water is
        flowing in the cell) divided by cell area. This would need to be updated
        to handle non-raster and/or non-uniform grids.

        Examples
        --------
        >>> from landlab import RasterModelGrid
        >>> from landlab.components import FlowAccumulator
        >>> grid = RasterModelGrid((3, 3), xy_spacing=1000.0)
        >>> elev = grid.add_zeros("topographic__elevation", at="node")
        >>> elev[3:] = 10.0
        >>> sed = grid.add_zeros("soil__depth", at="node")
        >>> sed[3:] = 100.0
        >>> fa = FlowAccumulator(grid)
        >>> fa.run_one_step()
        >>> eroder = GravelBedrockEroder(grid, abrasion_coefficient=0.0002)
        >>> eroder.calc_transport_rate()
        >>> eroder.calc_abrasion_rate()
        >>> int(eroder._abrasion[4] * 1e8)
        19
        """
        cores = self._grid.core_nodes
        self._abrasion[cores] = (
            self._abrasion_coef
            * 0.5
            * (self._sediment_outflux[cores] + self._sediment_influx[cores])
            * self._flow_link_length_over_cell_area
        )

    def calc_bedrock_abrasion_rate(self):
        """Update the rate of bedrock abrasion.

        Note: assumes _abrasion (of sediment) and _rock_exposure_fraction
        have already been updated. Like _abrasion, the rate is a length
        per time (equivalent to rate of lowering of the bedrock surface by
        abrasion). Result is stored in the field ``bedrock__abrasion_rate``.

        >>> import numpy as np
        >>> from landlab import RasterModelGrid
        >>> from landlab.components import FlowAccumulator
        >>> grid = RasterModelGrid((3, 4), xy_spacing=100.0)
        >>> elev = grid.add_zeros("topographic__elevation", at="node")
        >>> elev[:] = 0.01 * grid.x_of_node
        >>> sed = grid.add_zeros("soil__depth", at="node")
        >>> sed[:] = 1.0
        >>> fa = FlowAccumulator(grid)
        >>> fa.run_one_step()
        >>> eroder = GravelBedrockEroder(grid, abrasion_coefficient=1.0e-4)
        >>> eroder.calc_rock_exposure_fraction()
        >>> round(eroder._rock_exposure_fraction[6], 4)
        0.3679
        >>> eroder.calc_transport_rate()
        >>> np.round(eroder._sediment_outflux[5:7], 3)
        array([0.024, 0.012])
        >>> eroder.calc_abrasion_rate()
        >>> np.round(eroder._abrasion[5:7], 9)
        array([1.2e-08, 6.0e-09])
        >>> eroder.calc_bedrock_abrasion_rate()
        >>> np.round(eroder._rock_abrasion_rate[5:7], 10)
        array([4.4e-09, 2.2e-09])
        """
        self._rock_abrasion_rate = self._abrasion * self._rock_exposure_fraction

    def calc_bedrock_plucking_rate(self):
        """Update the rate of bedrock erosion by plucking.

        The rate is a volume per area per time [L/T], equivalent to the
        rate of lowering of the bedrock surface relative to the underlying
        material as a result of plucking. Result is stored in the field
        ``bedrock__plucking_rate``.

        Examples
        --------
        >>> import numpy as np
        >>> from landlab import RasterModelGrid
        >>> from landlab.components import FlowAccumulator
        >>> grid_res = 100.0
        >>> grid = RasterModelGrid((3, 3), xy_spacing=grid_res)
        >>> elev = grid.add_zeros("topographic__elevation", at="node")
        >>> elev[4] = 1.0
        >>> sed = grid.add_zeros("soil__depth", at="node")
        >>> fa = FlowAccumulator(grid)
        >>> fa.run_one_step()
        >>> eroder = GravelBedrockEroder(grid)
        >>> eroder.calc_rock_exposure_fraction()
        >>> eroder.calc_bedrock_plucking_rate()
        >>> predicted_plucking_rate = 1.0e-6 * 1.0e4 * 0.01 ** (7.0 / 6.0) / grid_res
        >>> round(predicted_plucking_rate, 9)  # Kp Q S^(7/6)
        4.64e-07
        >>> int(round(eroder._pluck_rate[4] * 1e9))
        464
        """
        cores = self._grid.core_nodes
        self._pluck_rate[cores] = (
            self._plucking_coef
            * self._intermittency_factor
            * self._discharge[cores]
            * self._slope[cores] ** _SEVEN_SIXTHS
            * self._rock_exposure_fraction[cores]
        ) * self._flow_link_length_over_cell_area

    def calc_sediment_influx(self):
        """Update the volume influx at each node.

        Result is stored in the field ``bedload_sediment__volume_influx``.
        """
        self._sediment_influx[:] = 0.0
        for c in self.grid.core_nodes:  # send sediment downstream
            r = self._receiver_node[c]
            self._sediment_influx[r] += self._sediment_outflux[c]

    def calc_sediment_rate_of_change(self):
        """Update the rate of thickness change of coarse sediment at each core node.

        Result is stored in the field ``sediment__rate_of_change``.

        Examples
        --------
        >>> import numpy as np
        >>> from landlab import RasterModelGrid
        >>> from landlab.components import FlowAccumulator
        >>> grid = RasterModelGrid((3, 4), xy_spacing=100.0)
        >>> elev = grid.add_zeros("topographic__elevation", at="node")
        >>> elev[:] = 0.01 * grid.x_of_node
        >>> sed = grid.add_zeros("soil__depth", at="node")
        >>> sed[:] = 100.0
        >>> grid.status_at_node[grid.perimeter_nodes] = grid.BC_NODE_IS_CLOSED
        >>> grid.status_at_node[4] = grid.BC_NODE_IS_FIXED_VALUE
        >>> fa = FlowAccumulator(grid)
        >>> fa.run_one_step()
        >>> eroder = GravelBedrockEroder(grid)
        >>> eroder.calc_transport_rate()
        >>> eroder.calc_sediment_influx()
        >>> eroder.calc_sediment_rate_of_change()
        >>> np.round(eroder._sediment_outflux[4:7], 3)
        array([0.   , 0.038, 0.019])
        >>> np.round(eroder._sediment_influx[4:7], 3)
        array([0.038, 0.019, 0.   ])
        >>> np.round(eroder._dHdt[5:7], 8)
        array([-2.93e-06, -2.93e-06])
        """
        cores = self.grid.core_nodes
        self._dHdt[cores] = self._porosity_factor * (
            (self._sediment_influx[cores] - self._sediment_outflux[cores])
            / self.grid.area_of_cell[self.grid.cell_at_node[cores]]
            + (self._pluck_rate[cores] * self._pluck_coarse_frac)
            - self._abrasion[cores]
        )

    def _update_slopes(self):
        """Update self._slope.

        Result is stored in field ``topographic__steepest_slope``.
        """
        dz = np.maximum(self._elev - self._elev[self._receiver_node], 0.0)
        if self._flow_length_is_variable:
            if self._grid_has_diagonals:
                link_len = self.grid.length_of_d8
            else:
                link_len = self.grid.length_of_link
            self._slope[self.grid.core_nodes] = (
                dz[self.grid.core_nodes] / link_len[self.grid.core_nodes]
            )
        else:
            self._slope[self.grid.core_nodes] = (
                dz[self.grid.core_nodes] / self.grid.spacing
            )

    def update_rates(self):
        """Update rate of sediment thickness change, and rate of bedrock lowering by abrasion
        and plucking.

        Combined rate of rock lowering relative to underlying material is stored in the field
        ``bedrock__lowering_rate``.

        Examples
        --------
        >>> import numpy as np
        >>> from landlab import RasterModelGrid
        >>> from landlab.components import FlowAccumulator
        >>> grid = RasterModelGrid((3, 4), xy_spacing=100.0)
        >>> elev = grid.add_zeros("topographic__elevation", at="node")
        >>> elev[:] = 0.01 * grid.x_of_node
        >>> sed = grid.add_zeros("soil__depth", at="node")
        >>> sed[:] = 1000.0
        >>> grid.status_at_node[grid.perimeter_nodes] = grid.BC_NODE_IS_CLOSED
        >>> grid.status_at_node[4] = grid.BC_NODE_IS_FIXED_VALUE
        >>> fa = FlowAccumulator(grid)
        >>> fa.run_one_step()
        >>> eroder = GravelBedrockEroder(grid)
        >>> eroder.run_one_step(1000.0)
        >>> np.round(elev[4:7], 4)
        array([0.    , 0.9971, 1.9971])
        """
        self._update_slopes()
        self.calc_rock_exposure_fraction()
        self.calc_transport_rate()
        self.calc_sediment_influx()

        if self._flow_length_is_variable:
            self._update_flow_link_length_over_cell_area()
        self.calc_bedrock_plucking_rate()
        if self._abrasion_coef > 0.0:
            self.calc_abrasion_rate()
            self.calc_bedrock_abrasion_rate()
        self.calc_sediment_rate_of_change()
        self._rock_lowering_rate = self._pluck_rate + self._rock_abrasion_rate

    def _update_rock_sed_and_elev(self, dt):
        """Update rock elevation, sediment thickness, and elevation
        using current rates of change extrapolated forward by time dt.
        """
        self._sed += self._dHdt * dt
        self._bedrock__elevation -= self._rock_lowering_rate * dt
        self._elev[:] = self._bedrock__elevation + self._sed

    def _estimate_max_time_step_size(self, upper_limit_dt=1.0e6):
        """
        Estimate the maximum possible time-step size that avoids
        flattening any streamwise slope or exhausting sediment.

        The ``upper_limit_dt`` parameter handles the special case of
        a nonexistent upper limit, which only occurs when there are
        no nodes at which either sediment or slope gradient is
        declining. Value is arbitrary as long as it is >= the user-provided
        global time-step size (in :meth:`~.run_one_step`).

        Parameters
        ----------
        dt : float (default 1.0e6)
            Maximum time step size
        """
        sed_is_declining = np.logical_and(self._dHdt < 0.0, self._sed > 0.0)
        if np.any(sed_is_declining):
            min_time_to_exhaust_sed = np.amin(
                -self._sed[sed_is_declining] / self._dHdt[sed_is_declining]
            )
        else:
            min_time_to_exhaust_sed = upper_limit_dt
        dzdt = self._dHdt - self._rock_lowering_rate
        rate_diff = dzdt[self._receiver_node] - dzdt
        height_above_rcvr = self._elev - self._elev[self._receiver_node]
        slope_is_declining = np.logical_and(rate_diff > 0.0, height_above_rcvr > 0.0)
        if np.any(slope_is_declining):
            min_time_to_flatten_slope = np.amin(
                height_above_rcvr[slope_is_declining] / rate_diff[slope_is_declining]
            )
        else:
            min_time_to_flatten_slope = upper_limit_dt
        return 0.5 * min(min_time_to_exhaust_sed, min_time_to_flatten_slope)

    def run_one_step(self, global_dt):
        """Advance solution by time interval global_dt, subdividing
        into sub-steps as needed."""
        time_remaining = global_dt
        while time_remaining > 0.0:
            self.update_rates()
            max_dt = self._estimate_max_time_step_size()
            this_dt = min(max_dt, time_remaining)
            this_dt = max(this_dt, _DT_MAX)
            self._update_rock_sed_and_elev(this_dt)
            time_remaining -= this_dt

#!/usr/bin/env python3
"""
Model bedrock incision and gravel transport and abrasion in a network of rivers.

@author: gtucker
"""

import numpy as np

from landlab import Component
from landlab import HexModelGrid
from landlab.grid.diagonals import DiagonalsMixIn

use_cfuncs = True
if use_cfuncs:
    from .cfuncs import _calc_sediment_influx
    from .cfuncs import _calc_sediment_rate_of_change
    from .cfuncs import _estimate_max_time_step_size_ext


_DT_MAX = 1.0e-2
_ONE_SIXTH = 1.0 / 6.0
_SEVEN_SIXTHS = 7.0 / 6.0
_D8_CHAN_LENGTH_FACTOR = 1.0  # 0.5 * (
#    1.0 + 2.0**0.5
# )  # D8 raster: average of straight and diagonal


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
    abrasion_coefficient : float (default 0.0 1/m) *DEPRECATED*
        Abrasion coefficient with units of inverse length
    sediment_porosity : float (default 0.35)
        Bulk porosity of bed sediment
    depth_decay_scale : float (default 1.0)
        Scale for depth decay in bedrock exposure function
    plucking_coefficient : float or (n_core_nodes,) array of float (default 1.0e-4 1/m)
        Rate coefficient for bedrock erosion by plucking
    number_of_sediment_classes : int (default 1)
        Number of sediment abradability classes
    init_thickness_per_class : float or (n_core_nodes,) array of float (default 1 / n-classes)
        Starting thickness for each sediment fraction
    abrasion_coefficients : iterable containing floats (default 0.0 1/m)
        Abrasion coefficients; should be same length as number of sed classes
    coarse_fractions_from_plucking : float or (n_core_nodes,) array of float (default 1.0)
        Fraction(s) of plucked material that becomes part of gravel sediment load
    rock_abrasion_index : int (default 0)
        If multiple classes, specifies which contains the abrasion
        coefficient for bedrock

    Notes
    -----
    The doctest below demonstrates approximate equilibrium between uplift, transport,
    and sediment abrasion in a case with effectively unlimited sediment. The
    analytical solution is:

    sediment input by uplift = sediment outflux + sediment loss to abrasion

    In math,

    U A = kq I Q S^(7/6) + 0.5 b Qs dx

    S = (U A / (kq I Q (1 + 0.5 b dx))) ^ 6/7 ~ 0.0342

    Elevation of the single core node = 0.0342 x 1,000 m ~ 34 m.
    However, because of VERY long time steps, the post-erosion elevation is 1 m
    lower, at 33 m (it will be uplifted by a meter at the start of each step).

    Examples
    --------
    >>> from landlab import RasterModelGrid
    >>> from landlab.components import FlowAccumulator
    >>> grid = RasterModelGrid((3, 3), xy_spacing=1000.0)
    >>> elev = grid.add_zeros("topographic__elevation", at="node")
    >>> sed = grid.add_zeros("soil__depth", at="node")
    >>> sed[4] = 1.0e6
    >>> grid.status_at_node[grid.perimeter_nodes] = grid.BC_NODE_IS_CLOSED
    >>> grid.status_at_node[5] = grid.BC_NODE_IS_FIXED_VALUE
    >>> fa = FlowAccumulator(grid, runoff_rate=10.0)
    >>> fa.run_one_step()
    >>> eroder = GravelBedrockEroder(
    ...     grid, sediment_porosity=0.0, abrasion_coefficients=[0.0005]
    ... )
    >>> rock_elev = grid.at_node["bedrock__elevation"]
    >>> for _ in range(150):
    ...     rock_elev[grid.core_nodes] += 1.0
    ...     elev[grid.core_nodes] += 1.0
    ...     fa.run_one_step()
    ...     eroder.run_one_step(10000.0)
    ...
    >>> int(elev[4])
    33
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
        sediment_porosity=0.35,
        depth_decay_scale=1.0,
        plucking_coefficient=1.0e-4,
        number_of_sediment_classes=1,
        init_thickness_per_class=None,
        abrasion_coefficients=0.0,
        coarse_fractions_from_plucking=1.0,
        rock_abrasion_index=0,
    ):
        """Initialize GravelBedrockEroder."""

        if init_thickness_per_class is None:
            init_thickness_per_class = 1.0 / number_of_sediment_classes
        abrasion_coefficients = np.broadcast_to(
            abrasion_coefficients, number_of_sediment_classes
        )
        init_thickness_per_class = np.broadcast_to(
            init_thickness_per_class, number_of_sediment_classes
        )
        coarse_fractions_from_plucking = np.broadcast_to(
            coarse_fractions_from_plucking, number_of_sediment_classes
        )

        super().__init__(grid)

        # Parameters
        self._trans_coef = transport_coefficient
        self._intermittency_factor = intermittency_factor
        self._porosity_factor = 1.0 / (1.0 - sediment_porosity)
        self._depth_decay_scale = depth_decay_scale
        if (
            isinstance(plucking_coefficient, np.ndarray)
            and len(plucking_coefficient) == self.grid.number_of_nodes
        ):
            plucking_coefficient = plucking_coefficient[self.grid.core_nodes]
        self._plucking_coef = plucking_coefficient
        self._rock_abrasion_index = rock_abrasion_index

        # Handle sediment classes, abrasion coefficients, and plucking fractions
        # if abrasion_coefficient is not None:
        #     if number_of_sediment_classes == 1:
        #         warn("Use abrasion_coefficients (plural)", DeprecationWarning)
        #         abrasion_coefficients = [abrasion_coefficient]
        #     else:
        #         warn(
        #             "Ignoring abrasion_coefficient; use abrasion_coefficients (plural)",
        #             DeprecationWarning,
        #         )
        # self._abrasion_coef = abrasion_coefficients[0]  # temporary for dev/test
        # if coarse_fraction_from_plucking is not None:
        #     if number_of_sediment_classes == 1:
        #         warn("Use coarse_fractions_from_plucking (plural)", DeprecationWarning)
        #         if (
        #             isinstance(coarse_fraction_from_plucking, np.ndarray)
        #             and len(coarse_fraction_from_plucking) == self.grid.number_of_nodes
        #         ):
        #             coarse_fraction_from_plucking = coarse_fraction_from_plucking[
        #                 self.grid.core_nodes
        #             ]
        #         coarse_fractions_from_plucking = [coarse_fraction_from_plucking]
        #     else:
        #         warn(
        #             "Ignoring abrasion_coefficient; use abrasion_coefficients (plural)",
        #             DeprecationWarning,
        #         )
        # self._pluck_coarse_frac = coarse_fractions_from_plucking[
        #    0
        # ]  # temporary for dev/test

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
        # self._abrasion = grid.at_node["bedload_sediment__rate_of_loss_to_abrasion"]
        self._rock_exposure_fraction = grid.at_node["bedrock__exposure_fraction"]
        self._rock_abrasion_rate = grid.at_node["bedrock__abrasion_rate"]
        self._pluck_rate = grid.at_node["bedrock__plucking_rate"]

        self._setup_length_of_flow_link()

        # Create 2d arrays
        self._thickness_by_class = np.zeros(
            (number_of_sediment_classes, grid.number_of_nodes)
        )
        for i in range(number_of_sediment_classes):
            self._thickness_by_class[i, :] = init_thickness_per_class[i] * self._sed

        self._sed_influxes = np.zeros(
            (number_of_sediment_classes, grid.number_of_nodes)
        )
        self._sed_outfluxes = np.zeros(
            (number_of_sediment_classes, grid.number_of_nodes)
        )
        self._sed_abr_rates = np.zeros(
            (number_of_sediment_classes, grid.number_of_nodes)
        )
        self._dHdt_by_class = np.zeros(
            (number_of_sediment_classes, grid.number_of_nodes)
        )
        if number_of_sediment_classes == 1:  # use view of array to save space
            self._sediment_outflux = self._sed_outfluxes[0]
            grid.at_node["bedload_sediment__volume_outflux"] = self._sediment_outflux
            self._sediment_influx = self._sed_influxes[0]
            grid.at_node["bedload_sediment__volume_influx"] = self._sediment_influx
        elif number_of_sediment_classes > 1:
            self._sediment_fraction = np.zeros(
                (number_of_sediment_classes, grid.number_of_nodes)
            )
        else:
            print("number_of_sediment_classes must be >0")
            raise ValueError

        self._num_sed_classes = number_of_sediment_classes
        self._abr_coefs = np.array(abrasion_coefficients)
        self._pluck_coarse_frac = coarse_fractions_from_plucking

    def _setup_length_of_flow_link(self):
        """Set up a float or array containing length of the flow link from
        each node, which is needed for the abrasion rate calculations.

        Note: if raster, assumes grid.dx == grid.dy
        """
        if isinstance(self.grid, HexModelGrid):
            self._flow_link_length_over_cell_area = (
                self.grid.spacing / self.grid.area_of_cell[0]
            )
            self._flow_length_is_variable = False
            self._grid_has_diagonals = False
        elif isinstance(self.grid, DiagonalsMixIn):
            self._flow_length_is_variable = False
            self._grid_has_diagonals = True
            self._flow_link_length_over_cell_area = (
                _D8_CHAN_LENGTH_FACTOR * self.grid.dx / self.grid.area_of_cell[0]
            )
            # self._update_flow_link_length_over_cell_area()
        else:
            self._flow_length_is_variable = True
            self._update_flow_link_length_over_cell_area()
            self._grid_has_diagonals = False

    def _update_flow_link_length_over_cell_area(self):
        """Update the ratio of the length of link along which water flows out of
        each node to the area of the node's cell."""
        self._flow_link_length_over_cell_area = (
            self.grid.length_of_link[self._receiver_link[self.grid.core_nodes]]
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

    def calc_sediment_fractions(self):
        """
        Calculate and store fraction of each sediment class in the sediment
        at each grid node.
        """
        self._sediment_fraction[:, :] = 0.0
        for i in range(self._num_sed_classes):
            self._sediment_fraction[i, :] = np.divide(
                self._thickness_by_class[i, :], self._sed, where=self._sed > 0.0
            )
        assert np.all(self._sed >= 0.0)  # temp test, to be removed

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
        array([ 0.,  1.])
        >>> sed[4] = 1.0  # exposure frac should be 1/e ~ 0.3679
        >>> sed[5] = 2.0  # exposure frac should be 1/e^2 ~ 0.1353
        >>> eroder.calc_rock_exposure_fraction()
        >>> np.round(eroder._rock_exposure_fraction[4:6], 4)
        array([ 0.3679,  0.1353])
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
        if self._num_sed_classes > 1:
            self.calc_sediment_fractions()
            for i in range(self._num_sed_classes):
                self._sed_outfluxes[i, :] = (
                    self._sediment_fraction[i, :] * self._sediment_outflux
                )

    def calc_abrasion_rate(self):
        """Update the volume rate of bedload loss to abrasion, per unit area.

        Here we use the average of incoming and outgoing sediment flux to
        calculate the loss rate to abrasion in each sediment class.
        The result is stored in self._sed_abr_rates.

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
        >>> eroder = GravelBedrockEroder(
        ...     grid,
        ...     number_of_sediment_classes=3,
        ...     abrasion_coefficients=[0.002, 0.0002, 0.00002],
        ... )
        >>> eroder.calc_transport_rate()
        >>> eroder.calc_abrasion_rate()
        >>> eroder._sed_abr_rates[:, 4]
        array([ 6.34350474e-07, 6.34350474e-08, 6.34350474e-09])
        """
        cores = self._grid.core_nodes
        for i in range(self._num_sed_classes):
            self._sed_abr_rates[i, cores] = (
                self._abr_coefs[i]
                * 0.5
                * (self._sed_outfluxes[i, cores] + self._sed_influxes[i, cores])
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
        >>> eroder = GravelBedrockEroder(grid, abrasion_coefficients=[1.0e-4])
        >>> eroder.calc_rock_exposure_fraction()
        >>> round(eroder._rock_exposure_fraction[6], 4)
        0.3679
        >>> eroder.calc_transport_rate()
        >>> np.round(eroder._sediment_outflux[5:7], 3)
        array([ 0.024,  0.012])
        >>> eroder.calc_abrasion_rate()
        >>> np.round(eroder._sed_abr_rates[0, 5:7], 9)
        array([  1.20000000e-08,   6.00000000e-09])
        >>> eroder.calc_bedrock_abrasion_rate()
        >>> np.round(eroder._rock_abrasion_rate[5:7], 10)
        array([  4.40000000e-09,   2.20000000e-09])
        """
        self._rock_abrasion_rate[:] = (
            self._sed_abr_rates[self._rock_abrasion_index]
            * self._rock_exposure_fraction
        )

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
        """
        Update the volume influx at each node for each sediment class.

        Results are stored:
        (1) Total flux in the field ``bedload_sediment__volume_influx``
        (2) Per size class in the _sed_influxes array
        """
        if use_cfuncs:
            _calc_sediment_influx(
                self._num_sed_classes,
                self.grid.number_of_core_nodes,
                self._sediment_influx,
                self._sed_influxes,
                self._sediment_outflux,
                self._sed_outfluxes,
                self.grid.core_nodes,
                self._receiver_node,
            )
        else:
            self._sediment_influx[:] = 0.0
            self._sed_influxes[:, :] = 0.0
            for c in self.grid.core_nodes:  # send sediment downstream
                r = self._receiver_node[c]
                if self._num_sed_classes > 1:
                    self._sediment_influx[r] += self._sediment_outflux[c]
                for i in range(self._num_sed_classes):
                    self._sed_influxes[i, r] += self._sed_outfluxes[i, c]

    def calc_sediment_rate_of_change(self):
        """
        Calculate and store time rate of change of sediment thickness
        by sediment class. Result stored in self._dHdt_by_class.

        The doctest below illustrates some of the steps in calculating
        the rate of change of sediment thickness. calc_transport_rate()
        sets the _sediment_outflux at each core node. The value of
        outflux can be calculated as follows:

        outflux = transport coefficient x intermittency factor
                  x discharge x slope^(7/6) x (1 - rock exposure fraction)

        Using default parameters, assuming a runoff rate of unity, and
        given the negligible bedrock exposure fraction,
        for the upstream grid cell this works out to:

        outflux = 0.041 x 0.01 x (100 x 100) x 0.01^(7/6) x 1
                ~ 0.19 m3/y

        For the downstream node, the discharge is doubled, so the
        flux is doubled.

        The call to calc_sediment_influx() then passes the outflux to
        each downstream neighbor (the "receiver" node), in the array
        _sediment_influx. For purposes of the test, we are *not*
        calculating sediment loss to abrasion, which therefore remains
        at its initial value of zero for each node.

        The call to calc_sediment_rate_of_change() adds up the influx
        and subtracts the outflux for each core node, dividing by its
        cell area and multiplying by the porosity factor to get a rate
        of thickness change. For both nodes, we have a net outflux of
        ~0.19 m3/y, a porosity factor of 1 / (1 - 0.35) ~ 1.54, and a
        grid cell area of 10,000 m2, so the rate of change should be
        ~2.93 x 10^-5 m/y. This is stored in _dHdt_by_class; because
        we have only on sediment class, these are the values that should
        be assigned to the single class at each of the two core nodes.

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
        >>> np.round(eroder._sed_outfluxes[0, 4:7], 3)
        array([ 0.   ,  0.038,  0.019])
        >>> np.round(eroder._sed_influxes[0, 4:7], 3)
        array([ 0.038,  0.019,  0.   ])
        >>> np.round(eroder._dHdt_by_class[0, 5:7], 8)
        array([ -2.93000000e-06,  -2.93000000e-06])
        """
        if use_cfuncs:
            _calc_sediment_rate_of_change(
                self._num_sed_classes,
                self.grid.number_of_core_nodes,
                self._porosity_factor,
                self.grid.area_of_cell[0],
                self.grid.core_nodes,
                self._pluck_coarse_frac,
                self._dHdt,
                self._pluck_rate,
                self._dHdt_by_class,
                self._sed_influxes,
                self._sed_outfluxes,
                self._sed_abr_rates,
            )
        else:
            cores = self.grid.core_nodes
            for i in range(self._num_sed_classes):
                self._dHdt_by_class[i, cores] = self._porosity_factor * (
                    (self._sed_influxes[i, cores] - self._sed_outfluxes[i, cores])
                    / self.grid.area_of_cell[self.grid.cell_at_node[cores]]
                    + (self._pluck_rate[cores] * self._pluck_coarse_frac[i])
                    - self._sed_abr_rates[i, cores]
                )
            self._dHdt[:] = np.sum(self._dHdt_by_class, axis=0)

    def _update_slopes(self):
        """Update self._slope.

        Result is stored in field ``topographic__steepest_slope``.
        """
        dz = np.maximum(self._elev - self._elev[self._receiver_node], 0.0)
        if self._flow_length_is_variable or self._grid_has_diagonals:
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

        The doctest below evaluates the code against the following
        calculation. We have two core nodes, each with a gradient of
        0.01 m/m and a surface area of 10,000 m2.

        The sediment cover is 0.69315 m, which happens to give a
        cover fraction of ~0.5.

        We run just one time step, in which the transport and abrasion
        rates apply to an initial slope of 0.01 and a discharge of
        10,000 and 20,000 m3/y at the upstream and downstream nodes,
        respectively.

        Transport rate, upstream node:

        coefficient x intermittency x discharge x slope^(7/6) x cover
        0.041       x 0.01          x 10,000 m3/y x 0.004642 x 0.5
        ~0.0095 m3/y

        Downstream node: ~0.019 m3/y (because of 2x discharge)

        Sediment influxes: at the open boundary node, should equal the
        outflux from the downstream core node (0.019 m3/y); at the
        downstream core node, should equal outflux from the upstream
        node (0.095 m3/y), and at the upstream node should be zero.

        With 50% bedrock exposure, the bedrock plucking rate for the
        upstream node should be:

        pluck coef x intermittency x discharge x slope^(7/6) x exposure
        / width of grid cell

        1e-4 x 0.01 x 10,000 x 0.01^(7/6) x 0.5 / 100 ~ 2.32e-7

        For the downstream node, it should be twice this value.

        With an abrasion coefficient of 0.001 1/m (fast!), the
        sediment abrasion rate at the upstream node should derive
        from the average of influx and outflux:

        0.5 x (0.0095 + 0) m3/y x 0.001 1/m  x 100 m ~ 0.000475 m3/y
        = 4.75e-8 m/y lowering equivalent

        For the downstream node,

        0.5 x (0.019 + 0.0095) m3/y x 0.001 1/m  x 100 m ~ 0.001425 m3/y
        = 1.425e-7 m/y lowering equivalent

        The sediment rate of change should be the sum of 3 terms:
        (1) the difference between influx and outflux, multiplied by
        the porosity factor 1 / (1 - phi) and divided by cell area.
        That equates to:
        -1 / (1 - 0.35) x 0.0095 / 10000 ~ -1.46 x 10^-6
        (2) the abrasion lowering rate above.
        (3) the plucking rate times the coarse fraction derived from
        plucking:
        + 2.32e-7 m/y x 0.25 ~ 5.8e-8 upstream, and 1.16e-7 downstream.
        Their sum is (approximately, subject to rounding errors):
        ~ -1.46e-6 - 1.43e-7 + 1.16e-7 = -1.487e-6 m/y downstream,
        ~ -1.46e-6 - 4.8e-8 + 5.8e-8 = -1.45e-6 m/y upstream.

        This should be stored in both _dHdt_by_class (the rate of
        thickening per sediment class, of which there's only one in
        this case) and _dHdt (which totals up all the classes).

        This rate is then extrapolated for one time step of 1000 y,
        so that the thickness is reduced by about 1.5 mm. The resulting
        thickness should be its original value minus this amount:
        0.69315 - 1.45e-3 = 0.6917 m, 0.69315 - 1.5e-3 = 0.69165 m

        Rock elevation should reduced by the loss of sediment plus the loss
        of rock. Rock loss is (plucking rate + rock abrasion rate) x time step.
        Upstream: (2 - 0.69315) - (2.32e-7 + 0.5 * 4.76e-8) m x 1000 y = 1.3065942 m
        Downstream: (1 - 0.69315) - (2.32e-7 + 0.5 * 4.76e-8) m x 1000 y = 0.3063145 m

        Finally, total elevation is the sum of the new sediment thickness
        and rock elevation (again, approximate, subject to rounding):
        Upstream: 0.3063 + 0.6917 = 0.9980
        Downstream: 1.3066 + 0.6917 = 1.9983

        Examples
        --------
        >>> import numpy as np
        >>> from landlab import RasterModelGrid
        >>> from landlab.components import FlowAccumulator
        >>> grid = RasterModelGrid((3, 4), xy_spacing=100.0)
        >>> elev = grid.add_zeros("topographic__elevation", at="node")
        >>> elev[:] = 0.01 * grid.x_of_node
        >>> sed = grid.add_zeros("soil__depth", at="node")
        >>> sed[:] = 0.69315
        >>> rock = grid.add_zeros("bedrock__elevation", at="node")
        >>> rock[:] = elev - sed
        >>> grid.status_at_node[grid.perimeter_nodes] = grid.BC_NODE_IS_CLOSED
        >>> grid.status_at_node[4] = grid.BC_NODE_IS_FIXED_VALUE
        >>> fa = FlowAccumulator(grid)
        >>> fa.run_one_step()
        >>> eroder = GravelBedrockEroder(
        ...     grid,
        ...     abrasion_coefficients=[0.001],
        ...     coarse_fractions_from_plucking=[0.25],
        ... )
        >>> eroder._bedrock__elevation[5:7]
        array([ 0.30685,  1.30685])
        >>> eroder.run_one_step(1000.0)
        >>> np.round(eroder._sed_outfluxes[0, 5:7], 3)
        array([ 0.019,  0.01 ])
        >>> np.round(eroder._sed_influxes[0, 4:7], 3)
        array([ 0.019,  0.01 ,  0.   ])
        >>> np.round(eroder._pluck_rate[5:7], 10)
        array([  4.64200000e-07,   2.32100000e-07])
        >>> np.round(eroder._sed_abr_rates[0, 5:7], 9)
        array([  1.43000000e-07,   4.80000000e-08])
        >>> np.round(eroder._dHdt_by_class[0, 5:7], 8)
        array([ -1.50000000e-06,  -1.45000000e-06])
        >>> np.round(eroder._dHdt[5:7], 8)
        array([ -1.50000000e-06,  -1.45000000e-06])
        >>> np.round(sed[5:7], 4)
        array([ 0.6916,  0.6917])
        >>> np.round(eroder._bedrock__elevation[5:7], 5)
        array([ 0.30631,  1.30659])
        >>> np.round(elev[4:7], 4)
        array([ 0.    ,  0.998 ,  1.9983])
        """
        self._update_slopes()
        self.calc_rock_exposure_fraction()
        self.calc_transport_rate()
        self.calc_sediment_influx()

        if self._flow_length_is_variable:
            self._update_flow_link_length_over_cell_area()
        self.calc_bedrock_plucking_rate()
        if np.amax(self._abr_coefs) > 0.0:
            self.calc_abrasion_rate()
            self.calc_bedrock_abrasion_rate()
        self.calc_sediment_rate_of_change()
        self._rock_lowering_rate[:] = self._pluck_rate + self._rock_abrasion_rate

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
        if use_cfuncs:
            min_dt = _estimate_max_time_step_size_ext(
                upper_limit_dt,
                self.grid.number_of_nodes,
                self._sed,
                self._elev,
                self._dHdt,
                self._dHdt - self._rock_lowering_rate,
                self._receiver_node,
            )
        else:
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
            slope_is_declining = np.logical_and(
                rate_diff > 0.0, height_above_rcvr > 0.0
            )
            if np.any(slope_is_declining):
                min_time_to_flatten_slope = np.amin(
                    height_above_rcvr[slope_is_declining]
                    / rate_diff[slope_is_declining]
                )
            else:
                min_time_to_flatten_slope = upper_limit_dt
            min_dt = 0.5 * min(min_time_to_exhaust_sed, min_time_to_flatten_slope)
        return min_dt

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

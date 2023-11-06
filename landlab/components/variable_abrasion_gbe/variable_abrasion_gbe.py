import numpy as np

from landlab.components import GravelBedrockEroder


class VariableAbrasionGBE(GravelBedrockEroder):
    """
    Landlab component that models fluvial landscape evolution based on
    Bedrock-Incising Gravel-Abrading Near-Threshold River (BIGANTR)
    theory, with the addition of allowing multiple sediment classes
    with independent abrasion coefficients.

    VariableAbrasionGBE is a subclass of GravelBedrockEroder.

    Parameters
    ----------
        grid : ModelGrid
        A Landlab model grid object
    intermittency_factor : float (default 0.01)
        Fraction of time that bankfull flow occurs
    transport_coefficient : float (default 0.041)
        Dimensionless transport efficiency factor (see Wickert & Schildgen 2019)
    number_of_sediment_classes : int (default 3)
        Number of sediment abradability classes
    abrasion_coefficients : iterable containing floats (default 0.0 1/m)
        Abrasion coefficients; should be same length as number of sed classes
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
    >>> grid.status_at_node[grid.perimeter_nodes] = grid.BC_NODE_IS_CLOSED
    >>> grid.status_at_node[5] = grid.BC_NODE_IS_FIXED_VALUE
    >>> elev[4] = 10.0
    >>> sed[4] = 300.0
    >>> fa = FlowAccumulator(grid)
    >>> fa.run_one_step()
    >>> vagbe = VariableAbrasionGBE(
    ...     grid,
    ...     sediment_porosity=0.0,
    ...     plucking_coefficient=0.0,
    ...     abrasion_coefficients=[0.002, 0.0002, 0.00002],
    ... )
    >>> vagbe._thickness_by_class[:, 4]
    array([ 100.,  100.,  100.])
    >>> vagbe.run_one_step(1.0)
    >>> grid.at_node["bedload_sediment__volume_outflux"][4]
    1.9030514217812389
    >>> grid.at_node["bedload_sediment__volume_influx"][5]
    1.9030514217812389
    >>> vagbe._sed_influxes[:, 5]
    array([ 0.63435047,  0.63435047,  0.63435047])
    >>> vagbe._sed_abr_rates[:, 4]
    array([ 6.34350474e-07, 6.34350474e-08, 6.34350474e-09])
    >>> vagbe._dHdt_by_class[:, 4]
    array([ -1.26870095e-06, -6.97785521e-07, -6.40693979e-07])
    """

    def __init__(
        self,
        grid,
        intermittency_factor=0.01,
        transport_coefficient=0.041,
        number_of_sediment_classes=3,
        init_thickness_per_class=None,
        abrasion_coefficients=0.0,
        sediment_porosity=0.35,
        depth_decay_scale=1.0,
        plucking_coefficient=1.0e-4,
        coarse_fractions_from_plucking=0.0,
    ):
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

        super().__init__(
            grid,
            intermittency_factor,
            transport_coefficient,
            max(abrasion_coefficients),
            sediment_porosity,
            depth_decay_scale,
            plucking_coefficient,
        )

        self._num_sed_classes = number_of_sediment_classes
        self._abr_coefs = abrasion_coefficients
        self._pluck_coarse_frac = coarse_fractions_from_plucking

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
        self._sediment_fraction = np.zeros(
            (number_of_sediment_classes, grid.number_of_nodes)
        )
        self._dHdt_by_class = np.zeros(
            (number_of_sediment_classes, grid.number_of_nodes)
        )

    def calc_sediment_fractions(self):
        """
        Calculate and store fraction of each sediment class in the sediment
        at each grid node.
        """
        for i in range(self._num_sed_classes):
            self._sediment_fraction[i, :] = self._thickness_by_class[i, :] / self._sed

    def calc_transport_rate(self):
        """

        Examples
        --------
        >>> from landlab import RasterModelGrid
        >>> from landlab.components import FlowAccumulator
        >>> grid = RasterModelGrid((3, 3), xy_spacing=1000.0)
        >>> elev = grid.add_zeros("topographic__elevation", at="node")
        >>> sed = grid.add_zeros("soil__depth", at="node")
        >>> grid.status_at_node[grid.perimeter_nodes] = grid.BC_NODE_IS_CLOSED
        >>> grid.status_at_node[5] = grid.BC_NODE_IS_FIXED_VALUE
        >>> elev[4] = 10.0
        >>> sed[4] = 300.0
        >>> fa = FlowAccumulator(grid)
        >>> fa.run_one_step()
        >>> vagbe = VariableAbrasionGBE(grid)
        >>> vagbe.run_one_step(1.0)
        >>> vagbe._sediment_fraction[:, 4]
        array([ 0.33333333,  0.33333333,  0.33333333])
        >>> vagbe._sed_outfluxes[:, 4]
        array([ 0.63435047,  0.63435047,  0.63435047])
        """
        super().calc_transport_rate()
        self.calc_sediment_fractions()
        for i in range(self._num_sed_classes):
            self._sed_outfluxes[i, :] = (
                self._sediment_fraction[i, :] * self._sediment_outflux
            )

    def calc_sediment_influx(self):
        """Update the volume influx at each node.

        Result is stored in the field ``bedload_sediment__volume_influx``.
        """
        self._sediment_influx[:] = 0.0
        self._sed_influxes[:, :] = 0.0
        for c in self.grid.core_nodes:  # send sediment downstream
            r = self._receiver_node[c]
            self._sediment_influx[r] += self._sediment_outflux[c]
            for i in range(self._num_sed_classes):
                self._sed_influxes[i, r] += self._sed_outfluxes[i, c]

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
        >>> eroder = VariableAbrasionGBE(
        ...     grid, abrasion_coefficients=[0.002, 0.0002, 0.00002]
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

    def calc_sediment_rate_of_change(self):
        cores = self.grid.core_nodes
        for i in range(self._num_sed_classes):
            self._dHdt_by_class[i, cores] = self._porosity_factor * (
                (self._sed_influxes[i, cores] - self._sed_outfluxes[i, cores])
                / self.grid.area_of_cell[self.grid.cell_at_node[cores]]
                + (self._pluck_rate[cores] * self._pluck_coarse_frac[i])
                - self._sed_abr_rates[i, cores]
            )

    def _update_rock_sed_and_elev(self, dt):
        """Update rock elevation, sediment thickness, and elevation
        using current rates of change extrapolated forward by time dt.
        """
        self._sed[:] = 0.0
        for i in range(self._num_sed_classes):
            self._thickness_by_class[i, :] += self._dHdt_by_class[i, :] * dt
            self._sed[:] += self._thickness_by_class[i, :]
        self._bedrock__elevation -= self._rock_lowering_rate * dt
        self._elev[:] = self._bedrock__elevation + self._sed

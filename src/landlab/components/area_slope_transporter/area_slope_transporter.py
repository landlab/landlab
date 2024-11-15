from landlab import Component


class AreaSlopeTransporter(Component):
    """Model drainage network evolution for a network of transport-limited
    rivers in which sediment transport rate is calculated as a power-law
    function of drainage area and local streamwise slope gradient.

    AreaSlopeTransporter is designed to operate together with a flow-routing
    component such as PriorityFloodFlowRouter, so that each grid node has
    a defined flow direction toward one of its neighbor nodes. Each core node
    is assumed to contain one outgoing fluvial channel, and (depending on
    the drainage structure) zero, one, or more incoming channels. These channels are
    treated as effectively sub-grid-scale features that are embedded in valleys
    that have a width of one grid cell. The rate of sediment transport out of
    a given node is calculated as a generic power function of drainage area,
    local slope, and a user-specified transport coefficient.
    Similar power-law formulations have been used, for example, by
    Willgoose et al. (1991a,b,c, and many papers following that use the
    SIBERIA model) and Howard (1994, in Water Resources Research).

    Parameters
    ----------
    grid : ModelGrid
        A Landlab model grid object
    transport_coefficient : float (default 0.0055)
        Dimensional transport efficiency factor
    area_exponent : float (default 1.4)
        Exponent on effective total discharge
    slope_exponent : float (default 2.1)
        Exponent on local streamwise slope gradient

    Examples
    --------
    >>> from landlab import RasterModelGrid
    >>> from landlab.components import FlowAccumulator
    >>> grid = RasterModelGrid((3, 3), xy_spacing=1000.0)
    >>> elev = grid.add_zeros("topographic__elevation", at="node")
    >>> grid.status_at_node[grid.perimeter_nodes] = grid.BC_NODE_IS_CLOSED
    >>> grid.status_at_node[5] = grid.BC_NODE_IS_FIXED_VALUE
    >>> fa = FlowAccumulator(grid)
    >>> fa.run_one_step()
    >>> transporter = AreaSlopeTransporter(grid)
    >>> for _ in range(200):
    ...     fa.run_one_step()
    ...     elev[grid.core_nodes] += 1.0
    ...     transporter.run_one_step(10000.0)
    ...
    >>> int(round(elev[4] * 100))
    1068
    """

    _name = "AreaSlopeTransporter"

    _unit_agnostic = True

    _info = {
        "sediment__volume_influx": {
            "dtype": float,
            "intent": "out",
            "optional": False,
            "units": "m**3/y",
            "mapping": "node",
            "doc": "Volumetric incoming streamwise sediment transport rate",
        },
        "sediment__volume_outflux": {
            "dtype": float,
            "intent": "out",
            "optional": False,
            "units": "m**3/y",
            "mapping": "node",
            "doc": "Volumetric outgoing streamwise sediment transport rate",
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
        "drainage_area": {
            "dtype": float,
            "intent": "in",
            "optional": False,
            "units": "m**2",
            "mapping": "node",
            "doc": "Upstream accumulated surface area contributing to the node's discharge",
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
        transport_coefficient=0.0055,
        area_exponent=1.4,
        slope_exponent=2.1,
    ):
        """Initialize AreaSlopeTransporter."""

        super().__init__(grid)

        # Parameters
        self._trans_coef = transport_coefficient
        self._area_exponent = area_exponent
        self._slope_exponent = slope_exponent

        # Fields and arrays
        self._elev = grid.at_node["topographic__elevation"]
        self._area = grid.at_node["drainage_area"]
        self._slope = grid.at_node["topographic__steepest_slope"]
        self._receiver_node = grid.at_node["flow__receiver_node"]
        self._receiver_link = grid.at_node["flow__link_to_receiver_node"]
        super().initialize_output_fields()
        self._sediment_influx = grid.at_node["sediment__volume_influx"]
        self._sediment_outflux = grid.at_node["sediment__volume_outflux"]
        self._dzdt = grid.at_node["sediment__rate_of_change"]

    def calc_transport_capacity(self):
        """Calculate and return bed-load transport capacity.

        Calculation uses power-law approach, and provides
        volume per time rate.

        Examples
        --------
        >>> from landlab import RasterModelGrid
        >>> from landlab.components import FlowAccumulator
        >>> grid = RasterModelGrid((3, 3), xy_spacing=100.0)
        >>> elev = grid.add_zeros("topographic__elevation", at="node")
        >>> elev[3:] = 1.0
        >>> fa = FlowAccumulator(grid)
        >>> fa.run_one_step()
        >>> transporter = AreaSlopeTransporter(grid)
        >>> transporter.calc_transport_capacity()
        >>> int(transporter._sediment_outflux[4] * 1000)
        138
        """
        self._sediment_outflux[:] = (
            self._trans_coef
            * self._area**self._area_exponent
            * self._slope**self._slope_exponent
        )

    def calc_sediment_rate_of_change(self):
        """Update the rate of thickness change of sediment at each core node.

        Examples
        --------
        >>> import numpy as np
        >>> from landlab import RasterModelGrid
        >>> from landlab.components import FlowAccumulator
        >>> grid = RasterModelGrid((3, 4), xy_spacing=100.0)
        >>> elev = grid.add_zeros("topographic__elevation", at="node")
        >>> elev[:] = 0.01 * grid.x_of_node
        >>> grid.status_at_node[grid.perimeter_nodes] = grid.BC_NODE_IS_CLOSED
        >>> grid.status_at_node[4] = grid.BC_NODE_IS_FIXED_VALUE
        >>> fa = FlowAccumulator(grid)
        >>> fa.run_one_step()
        >>> transporter = AreaSlopeTransporter(grid)
        >>> transporter.calc_sediment_rate_of_change()
        >>> np.round(transporter._sediment_outflux[4:7], 3)
        array([0.   , 0.365, 0.138])
        >>> np.round(transporter._sediment_influx[4:7], 3)
        array([0.365, 0.138, 0.   ])
        >>> np.round(transporter._dzdt[5:7], 8)
        array([-2.264e-05, -1.382e-05])
        """
        self.calc_transport_capacity()
        cores = self.grid.core_nodes
        self._sediment_influx[:] = 0.0
        for c in cores:  # send sediment downstream
            r = self._receiver_node[c]
            self._sediment_influx[r] += self._sediment_outflux[c]
        self._dzdt[cores] = (
            self._sediment_influx[cores] - self._sediment_outflux[cores]
        ) / self.grid.area_of_cell[self.grid.cell_at_node[cores]]

    def run_one_step(self, dt):
        """Advance solution by time interval dt.

        Examples
        --------
        >>> import numpy as np
        >>> from landlab import RasterModelGrid
        >>> from landlab.components import FlowAccumulator
        >>> grid = RasterModelGrid((3, 4), xy_spacing=100.0)
        >>> elev = grid.add_zeros("topographic__elevation", at="node")
        >>> elev[:] = 0.01 * grid.x_of_node
        >>> grid.status_at_node[grid.perimeter_nodes] = grid.BC_NODE_IS_CLOSED
        >>> grid.status_at_node[4] = grid.BC_NODE_IS_FIXED_VALUE
        >>> fa = FlowAccumulator(grid)
        >>> fa.run_one_step()
        >>> transporter = AreaSlopeTransporter(grid)
        >>> transporter.run_one_step(10000.0)
        >>> np.round(elev[4:7], 4)
        array([0.    , 0.7736, 1.8618])
        """
        self.calc_sediment_rate_of_change()
        self._elev += self._dzdt * dt

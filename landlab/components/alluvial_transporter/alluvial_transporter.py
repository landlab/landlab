import numpy as np
from landlab import Component


class GravelRiverTransporter(Component):

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
            "units": "m**2/y",
            "mapping": "node",
            "doc": "Volumetric incoming streamwise bedload sediment transport rate",
        },
        "bedload_sediment__volume_outflux": {
            "dtype": float,
            "intent": "out",
            "optional": False,
            "units": "m**2/y",
            "mapping": "node",
            "doc": "Volumetric outgoing streamwise bedload sediment transport rate",
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
    ):
        """Initialize AlluvialTransporter1."""

        super().__init__(grid)

        # Parameters
        self._trans_coef = transport_coefficient
        self._intermittency_factor = intermittency_factor
        self._abrasion_coef = abrasion_coefficient

        # Fields and arrays
        self._elev = grid.at_node["topographic__elevation"]
        self._discharge = grid.at_node["surface_water__discharge"]
        self._slope = grid.at_node["topographic__steepest_slope"]
        self._receiver_node = grid.at_node["flow__receiver_node"]
        self._receiver_link = grid.at_node["flow__link_to_receiver_node"]
        super().initialize_output_fields()
        self._sediment_influx = grid.at_node["bedload_sediment__volume_influx"]
        self._sediment_outflux = grid.at_node["bedload_sediment__volume_outflux"]
        self._dzdt = grid.at_node["sediment__rate_of_change"]
        self._abrasion = grid.at_node["bedload_sediment__rate_of_loss_to_abrasion"]

        # Constants
        self._SEVEN_SIXTHS = 7.0 / 6.0

    def calc_transport_capacity(self):
        """Calculate and return bed-load transport capacity.

        Calculation uses Wickert-Schildgen approach, and provides
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
        >>> transporter = AlluvialTransporter1(grid)
        >>> transporter.calc_transport_capacity()
        >>> round(transporter._sediment_outflux[4], 4)
        0.019
        """
        self._sediment_outflux[:] = (
            self._trans_coef
            * self._intermittency_factor
            * self._discharge
            * self._slope**self._SEVEN_SIXTHS
        )

    def calc_abrasion_rate(self):
        """Update the rate of bedload loss to abrasion, per unit area."""
        cores = self._grid.core_nodes
        self._abrasion[cores] = (
            self._abrasion_coef
            * self._sediment_outflux[cores]
            / self._grid.length_of_link[self._receiver_link[cores]]
        )

    def calc_sediment_rate_of_change(self):
        """Update the rate of thickness change of coarse sediment at each core node.

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
        >>> transporter = AlluvialTransporter1(grid)
        >>> transporter.calc_sediment_rate_of_change()
        >>> np.round(transporter._sediment_outflux[4:7], 3)
        array([ 0.   ,  0.038,  0.019])
        >>> np.round(transporter._sediment_influx[4:7], 3)
        array([ 0.038,  0.019,  0.   ])
        >>> np.round(transporter._dzdt[5:7], 8)
        array([ -1.90000000e-06,  -1.90000000e-06])
        """
        self.calc_transport_capacity()
        cores = self.grid.core_nodes
        self._sediment_influx[self._receiver_node[cores]] = self._sediment_outflux[
            cores
        ]
        self._dzdt[cores] = (
            (self._sediment_influx[cores] - self._sediment_outflux[cores])
            / self.grid.area_of_cell[self.grid.cell_at_node[cores]]
        ) - self._abrasion

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
        >>> transporter = AlluvialTransporter1(grid)
        >>> transporter.run_one_step(1000.0)
        >>> np.round(elev[4:7], 4)
        array([ 0.    ,  0.9981,  1.9981])
        """
        self.calc_sediment_rate_of_change()
        self._elev += self._dzdt * dt

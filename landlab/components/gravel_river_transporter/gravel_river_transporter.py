import numpy as np
from scipy.sparse.linalg import spsolve

from landlab import Component
from landlab import HexModelGrid
from landlab.grid.diagonals import DiagonalsMixIn


def make_empty_matrix_and_rhs(grid):
    from scipy.sparse import csc_matrix

    mat = csc_matrix((grid.number_of_core_nodes, grid.number_of_core_nodes))
    rhs = np.zeros(grid.number_of_core_nodes)
    return mat, rhs


def zero_out_matrix(grid, mat, rcvr, mat_id):
    for i in grid.core_nodes:
        j = mat_id[i]
        mat[j, j] = 0.0
        r = rcvr[i]
        if grid.status_at_node[r] == grid.BC_NODE_IS_CORE:
            k = mat_id[r]
            mat[j, k] = 0.0
            mat[k, k] = 0.0
            mat[k, j] = 0.0


class GravelRiverTransporter(Component):
    """Model drainage network evolution for a network of transport-limited
    gravel-bed rivers with downstream abrasion.

    GravelRiverTransporter is designed to operate together with a flow-routing
    component such as PriorityFloodFlowRouter, so that each grid node has
    a defined flow direction toward one of its neighbor nodes. Each core node
    is assumed to contain one outgoing fluvial channel, and (depending on
    the drainage structure) zero, one, or more incoming channels. These channels are
    treated as effectively sub-grid-scale features that are embedded in valleys
    that have a width of one grid cell. The rate of gravel transport out of
    a given node is calculated as the product of bankfull discharge, channel
    gradient (to the 7/6 power), a dimensionless transport coefficient, and
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
    distance (for example, if a grid node is carrying a gravel load Qs to a
    neighboring node dx meters downstream, the rate of gravel loss in volume per
    time per area at the node will be beta Qs dx, where beta is the abrasion
    coefficient). Sediment mass conservation is calculated across each entire
    grid cell. For example, if a cell has surface area A, a total volume influx
    Qin, and downstream transport rate Qs, the resulting rate of change of
    elevation will be (Qin - Qs / (A (1 - phi)), where phi is porosity.

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
    solver : string (default "explicit")
        Solver type (currently only "explicit" is tested and operational)

    Examples
    --------
    >>> from landlab import RasterModelGrid
    >>> from landlab.components import FlowAccumulator
    >>> grid = RasterModelGrid((3, 3), xy_spacing=1000.0)
    >>> elev = grid.add_zeros("topographic__elevation", at="node")
    >>> grid.status_at_node[grid.perimeter_nodes] = grid.BC_NODE_IS_CLOSED
    >>> grid.status_at_node[5] = grid.BC_NODE_IS_FIXED_VALUE
    >>> fa = FlowAccumulator(grid, runoff_rate=10.0)
    >>> fa.run_one_step()
    >>> transporter = GravelRiverTransporter(grid, abrasion_coefficient=0.0005)
    >>> for _ in range(200):
    ...     fa.run_one_step()
    ...     elev[grid.core_nodes] += 1.0
    ...     transporter.run_one_step(10000.0)
    ...
    >>> int(elev[4] * 100)
    2366
    """

    _ONE_SIXTH = 1.0 / 6.0

    _name = "GravelRiverTransporter"

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
        sediment_porosity=0.35,
        solver="explicit",
    ):
        """Initialize GravelRiverTransporter."""

        super().__init__(grid)

        # Parameters
        self._trans_coef = transport_coefficient
        self._intermittency_factor = intermittency_factor
        self._abrasion_coef = abrasion_coefficient
        self._porosity_factor = 1.0 / (1.0 - sediment_porosity)

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

        # Solver type
        if solver == "explicit":
            self.run_one_step = self.run_one_step_simple_explicit
        elif solver == "matrix":
            import warnings

            from landlab.utils.matrix import get_core_node_at_node

            warnings.warn(
                "Matrix-based solver is experimental & not fully tested", stacklevel=2
            )
            self.run_one_step = self.run_one_step_matrix_inversion
            self._mat, self._rhs = make_empty_matrix_and_rhs(grid)
            self._mat_id = np.zeros(grid.number_of_nodes, dtype=int)
            self._mat_id = get_core_node_at_node(grid)
        else:
            raise ValueError("Solver type not recognized")

        self._setup_length_of_flow_link()

    def _setup_length_of_flow_link(self):
        """Set up a float or array containing length of the flow link from each node,
        which is needed for the abrasion rate calculations.
        """
        if isinstance(self.grid, HexModelGrid):
            self._flow_link_length_over_cell_area = (
                self.grid.spacing / self.grid.area_of_cell[0]
            )
            self._flow_length_is_variable = False
        elif isinstance(self.grid, DiagonalsMixIn):
            self._flow_length_is_variable = True
            self._grid_has_diagonals = True
        else:
            self._flow_length_is_variable = True
            self._grid_has_diagonals = False

    def _update_flow_link_length_over_cell_area(self):
        """Update the ratio of the length of link along which water flows out of
        each node to the area of the node's cell."""
        if self._grid_has_diagonals:
            self._flow_link_length_over_cell_area = (
                self.grid.length_of_d8[self._receiver_link[self.grid.core_nodes]]
                / self.grid.area_of_cell[self.grid.cell_at_node[self.grid.core_nodes]]
            )
        else:
            self._flow_link_length_over_cell_area = (
                self.grid.length_of_link[self._receiver_link[self.grid.core_nodes]]
                / self.grid.area_of_cell[self.grid.cell_at_node[self.grid.core_nodes]]
            )

    def calc_implied_depth(self, grain_diameter=0.01):
        """Utility function that calculates and returns water depth implied by
        slope and grain diameter, using Wickert & Schildgen (2019) equation 8.

        The equation is

            h = ((rho_s - rho / rho)) (1 + epsilon) tau_c* (D / S)

        where the factors on the right are sediment and water density, excess
        shear-stress factor, critical Shields stress, grain diameter, and slope
        gradient. Here the prefactor on D/S assumes sediment density of 2650 kg/m3,
        water density of 1000 kg/m3, shear-stress factor of 0.2, and critical
        Shields stress of 0.0495, giving a value of 0.09801.

        Examples
        --------
        >>> from landlab import RasterModelGrid
        >>> from landlab.components import FlowAccumulator
        >>> grid = RasterModelGrid((3, 3), xy_spacing=1000.0)
        >>> elev = grid.add_zeros("topographic__elevation", at="node")
        >>> elev[3:] = 10.0
        >>> fa = FlowAccumulator(grid)
        >>> fa.run_one_step()
        >>> transporter = GravelRiverTransporter(grid)
        >>> depth = transporter.calc_implied_depth()
        >>> int(depth[4] * 1000)
        98
        """
        DEPTH_FACTOR = 0.09801
        depth = np.zeros(self._grid.number_of_nodes)
        nonzero_slope = self._slope > 0.0
        depth[nonzero_slope] = (
            DEPTH_FACTOR * grain_diameter / self._slope[nonzero_slope]
        )
        return depth

    def calc_implied_width(self, grain_diameter=0.01, time_unit="y"):
        """Utility function that calculates and returns channel width implied by
        discharge, slope, and grain diameter, using Wickert & Schildgen (2019)
        equation 16.

        The equation is

            b = kb Q S**(7/6) / D**(3/2)

        where the dimensional prefactor, which includes sediment and water
        density, gravitational acceleration, critical Shields stress, and the
        transport factor epsilon, is

            kb = 0.17 g**(-1/2) (((rho_s - rho) / rho) (1 + eps) tau_c*)**(-5/3)

        Using g = 9.8 m/s2, rho_s = 2650 (quartz), rho = 1000 kg/m3, eps = 0.2,
        and tau_c* = 0.0495, kb ~ 2.61 s/m**(1/2). Converting to years,
        kb = 8.26e-8.

        Examples
        --------
        >>> from landlab import RasterModelGrid
        >>> from landlab.components import FlowAccumulator
        >>> grid = RasterModelGrid((3, 3), xy_spacing=10000.0)
        >>> elev = grid.add_zeros("topographic__elevation", at="node")
        >>> elev[3:] = 100.0
        >>> fa = FlowAccumulator(grid)
        >>> fa.run_one_step()
        >>> transporter = GravelRiverTransporter(grid)
        >>> width = transporter.calc_implied_width()
        >>> int(width[4] * 100)
        3833
        >>> grid.at_node["surface_water__discharge"] *= 1.0 / (3600 * 24 * 365.25)
        >>> width = transporter.calc_implied_width(time_unit="s")
        >>> int(width[4] * 100)
        3838
        """
        if time_unit[0] == "y":
            width_fac = 8.26e-8
        else:
            width_fac = 2.61  # assume seconds if not years
        width = np.zeros(self._grid.number_of_nodes)
        width = (
            width_fac
            * self._discharge
            * self._slope ** (7.0 / 6.0)
            / (grain_diameter**1.5)
        )
        return width

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
        >>> transporter = GravelRiverTransporter(grid)
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
        """Update the rate of bedload loss to abrasion, per unit area.

        Here we use the average of incoming and outgoing sediment flux to
        calculate the loss rate to abrasion.

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
        >>> fa = FlowAccumulator(grid)
        >>> fa.run_one_step()
        >>> transporter = GravelRiverTransporter(grid, abrasion_coefficient=0.0002)
        >>> transporter.calc_transport_capacity()
        >>> transporter.calc_abrasion_rate()
        >>> int(transporter._abrasion[4] * 1e8)
        19
        """
        cores = self._grid.core_nodes
        if self._flow_length_is_variable:
            self._update_flow_link_length_over_cell_area()
        self._abrasion[cores] = (
            self._abrasion_coef
            * 0.5
            * (self._sediment_outflux[cores] + self._sediment_influx[cores])
            * self._flow_link_length_over_cell_area
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
        >>> transporter = GravelRiverTransporter(grid)
        >>> transporter.calc_sediment_rate_of_change()
        >>> np.round(transporter._sediment_outflux[4:7], 3)
        array([0.   , 0.038, 0.019])
        >>> np.round(transporter._sediment_influx[4:7], 3)
        array([0.038, 0.019, 0.   ])
        >>> np.round(transporter._dzdt[5:7], 8)
        array([-2.93e-06, -2.93e-06])
        """
        self.calc_transport_capacity()
        if self._abrasion_coef > 0.0:
            self.calc_abrasion_rate()
        cores = self.grid.core_nodes
        self._sediment_influx[:] = 0.0
        for c in cores:  # send sediment downstream
            r = self._receiver_node[c]
            self._sediment_influx[r] += self._sediment_outflux[c]
        self._dzdt[cores] = self._porosity_factor * (
            (self._sediment_influx[cores] - self._sediment_outflux[cores])
            / self.grid.area_of_cell[self.grid.cell_at_node[cores]]
            - self._abrasion[cores]
        )

    def run_one_step_simple_explicit(self, dt):
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
        >>> transporter = GravelRiverTransporter(grid, solver="explicit")
        >>> transporter.run_one_step(1000.0)
        >>> np.round(elev[4:7], 4)
        array([0.    , 0.9971, 1.9971])
        """
        self.calc_sediment_rate_of_change()
        self._elev += self._dzdt * dt

    def _fill_matrix_and_rhs(self, dt):
        """Fill out entries in a sparse matrix and corresponding right-hand side
        vector.

        Examples
        --------
        >>> from landlab import RasterModelGrid
        >>> from landlab.components import FlowAccumulator
        >>> grid = RasterModelGrid((3, 4), xy_spacing=100.0)
        >>> elev = grid.add_zeros("topographic__elevation", at="node")
        >>> elev[:] = 0.01 * grid.x_of_node
        >>> grid.status_at_node[grid.perimeter_nodes] = grid.BC_NODE_IS_CLOSED
        >>> grid.status_at_node[4] = grid.BC_NODE_IS_FIXED_VALUE
        >>> fa = FlowAccumulator(grid)
        >>> transporter = GravelRiverTransporter(grid, solver="matrix")
        >>> transporter._mat.toarray()
        array([[0., 0.],
               [0., 0.]])
        >>> fa.run_one_step()
        >>> transporter._receiver_node[5:7]
        array([4, 5])
        >>> transporter._fill_matrix_and_rhs(1000.0)
        >>> transporter._rhs
        array([1., 2.])
        """
        prefac = (
            self._trans_coef * self._intermittency_factor * self._porosity_factor * dt
        ) / self.grid.dx**2
        a = prefac * (1.0 / self.grid.dx + self._abrasion_coef / 2)
        b = prefac * (1.0 / self.grid.dx - self._abrasion_coef / 2)
        f = self._discharge * (self._slope ** (self._ONE_SIXTH))

        zero_out_matrix(self.grid, self._mat, self._receiver_node, self._mat_id)
        for i in self.grid.core_nodes:
            j = self._mat_id[i]
            self._rhs[j] = self._elev[i]
            self._mat[j, j] += 1 + a * f[i]
            r = self._receiver_node[i]
            if self.grid.status_at_node[r] == self.grid.BC_NODE_IS_CORE:
                k = self._mat_id[r]
                self._mat[j, k] -= a * f[i]
                self._mat[k, k] += b * f[i]
                self._mat[k, j] -= b * f[i]
            else:
                self._rhs[j] += a * f[i] * self._elev[r]

    def run_one_step_matrix_inversion(self, dt):
        """Advance solution by time interval dt.

        WARNING: EXPERIMENTAL AND NOT FULLY TESTED - USE AT OWN RISK!

        Notes
        -----
        Does not update abrasion rate or sediment outflux fields.

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
        >>> transporter = GravelRiverTransporter(grid, solver="matrix")
        >>> transporter.run_one_step == transporter.run_one_step_matrix_inversion
        True
        >>> transporter.run_one_step(1000.0)
        """
        self._fill_matrix_and_rhs(dt)
        self._elev[self.grid.core_nodes] = spsolve(self._mat, self._rhs)

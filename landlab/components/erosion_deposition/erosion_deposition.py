import numpy as np

from landlab.utils.return_array import return_array_at_node

from ..depression_finder.lake_mapper import _FLOODED
from .cfuncs import calculate_qs_in
from .generalized_erosion_deposition import _GeneralizedErosionDeposition

ROOT2 = np.sqrt(2.0)  # syntactic sugar for precalculated square root of 2
TIME_STEP_FACTOR = 0.5  # factor used in simple subdivision solver
DEFAULT_MINIMUM_TIME_STEP = 0.001  # default minimum time step duration


class ErosionDeposition(_GeneralizedErosionDeposition):
    r"""
    Erosion-Deposition model in the style of Davy and Lague (2009). It uses a
    mass balance approach across the total sediment mass both in the bed and
    in transport coupled with explicit representation of the sediment
    transport lengthscale (the "xi-q" model) to derive a range of erosional
    and depositional responses in river channels.

    This implementation is close to the Davy & Lague scheme, with a few
    deviations:

        - Sediment porosity is handled explicitly in this implementation.

        - A fraction of the eroded sediment is permitted to enter the wash load,
          and lost to the mass balance (`F_f`).

        - Here an incision threshold :math:`\omega` is permitted, where it was not by Davy &
          Lague. It is implemented with an exponentially smoothed form to prevent
          discontinuities in the parameter space. See the
          :py:class:`~landlab.components.StreamPowerSmoothThresholdEroder`
          for more documentation.

        - This component uses an "effective" settling velocity, v_s, as one of its
          inputs. This parameter is simply equal to Davy & Lague's `d_star * V`
          dimensionless number.

    Erosion of the bed follows a stream power formulation, i.e.,

    .. math:

        E = K * q ** m_{sp} * S ** {n_sp} - \omega

    Note that the transition between transport-limited and detachment-limited
    behavior is controlled by the dimensionless ratio (v_s/r) where r is the
    runoff ratio (Q=Ar). r can be changed in the flow accumulation component
    but is not changed within ErosionDeposition. Because the runoff ratio r
    is not changed within the ErosionDeposition component,  v_s becomes the
    parameter that fundamentally controls response style. Very small v_s will
    lead to a detachment-limited response style, very large v_s will lead to a
    transport-limited response style. v_s == 1 means equal contributions from
    transport and erosion, and a hybrid response as described by Davy & Lague.

    Component written by C. Shobe, K. Barnhart, and G. Tucker.
    """

    _name = "ErosionDeposition"

    _info = {
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
        "sediment__flux": {
            "dtype": float,
            "intent": "out",
            "optional": False,
            "units": "m3/s",
            "mapping": "node",
            "doc": "Sediment flux (volume per unit time of sediment entering each node)",
        },
        "surface_water__discharge": {
            "dtype": float,
            "intent": "in",
            "optional": False,
            "units": "m**2/s",
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
        K=0.002,
        phi=0.3,
        v_s=1.0,
        m_sp=0.5,
        n_sp=1.0,
        sp_crit=0.0,
        F_f=0.0,
        discharge_field="surface_water__discharge",
        solver="basic",
        erode_flooded_nodes=True,
        dt_min=DEFAULT_MINIMUM_TIME_STEP,
    ):
        """Initialize the ErosionDeposition model.

        Parameters
        ----------
        grid : ModelGrid
            Landlab ModelGrid object
        K : float, field name, or array
            Erodibility for substrate (units vary).
        phi : float
            Sediment porosity [-].
        v_s : float
            Effective settling velocity for chosen grain size metric [L/T].
        m_sp : float
            Discharge exponent (units vary)
        n_sp : float
            Slope exponent (units vary)
        sp_crit : float, field name, or array
            Critical stream power to erode substrate [E/(TL^2)]
        F_f : float
            Fraction of eroded material that turns into "fines" that do not
            contribute to (coarse) sediment load. Defaults to zero.
        discharge_field : float, field name, or array
            Discharge [L^2/T]. The default is to use the grid field
            'surface_water__discharge', which is simply drainage area
            multiplied by the default rainfall rate (1 m/yr). To use custom
            spatially/temporally varying rainfall, use 'water__unit_flux_in'
            to specify water input to the FlowAccumulator.
        erode_flooded_nodes : bool (optional)
            Whether erosion occurs in flooded nodes identified by a
            depression/lake mapper (e.g., DepressionFinderAndRouter). When set
            to false, the field *flood_status_code* must be present on the grid
            (this is created by the DepressionFinderAndRouter). Default True.
        solver : string
            Solver to use. Options at present include:
                (1) 'basic' (default): explicit forward-time extrapolation.
                    Simple but will become unstable if time step is too large.
                (2) 'adaptive': adaptive time-step solver that estimates a
                    stable step size based on the shortest time to "flattening"
                    among all upstream-downstream node pairs.

        Examples
        ---------
        >>> import numpy as np
        >>> from landlab import RasterModelGrid
        >>> from landlab.components import FlowAccumulator
        >>> from landlab.components import DepressionFinderAndRouter
        >>> from landlab.components import ErosionDeposition
        >>> from landlab.components import FastscapeEroder
        >>> np.random.seed(seed = 5000)

        Define grid and initial topography:
            -5x5 grid with baselevel in the lower left corner
            -all other boundary nodes closed
            -Initial topography is plane tilted up to the upper right + noise

        >>> nr = 5
        >>> nc = 5
        >>> dx = 10
        >>> mg = RasterModelGrid((nr, nc), xy_spacing=10.0)
        >>> _ = mg.add_zeros('node', 'topographic__elevation')
        >>> mg['node']['topographic__elevation'] += (mg.node_y/10 +
        ...        mg.node_x/10 + np.random.rand(len(mg.node_y)) / 10)
        >>> mg.set_closed_boundaries_at_grid_edges(bottom_is_closed=True,
        ...                                               left_is_closed=True,
        ...                                               right_is_closed=True,
        ...                                               top_is_closed=True)
        >>> mg.set_watershed_boundary_condition_outlet_id(0,\
                mg['node']['topographic__elevation'], -9999.)
        >>> fsc_dt = 100.
        >>> ed_dt = 1.

        Check initial topography

        >>> mg.at_node['topographic__elevation'] # doctest: +NORMALIZE_WHITESPACE
        array([ 0.02290479,  1.03606698,  2.0727653 ,  3.01126678,  4.06077707,
            1.08157495,  2.09812694,  3.00637448,  4.07999597,  5.00969486,
            2.04008677,  3.06621577,  4.09655859,  5.04809001,  6.02641123,
            3.05874171,  4.00585786,  5.0595697 ,  6.04425233,  7.05334077,
            4.05922478,  5.0409473 ,  6.07035008,  7.0038935 ,  8.01034357])

        Instantiate Fastscape eroder, flow router, and depression finder

        >>> fr = FlowAccumulator(mg, flow_director='D8')
        >>> df = DepressionFinderAndRouter(mg)
        >>> fsc = FastscapeEroder(
        ...     mg,
        ...     K_sp=.001,
        ...     m_sp=.5,
        ...     n_sp=1,
        ...     erode_flooded_nodes=False)

        Burn in an initial drainage network using the Fastscape eroder:

        >>> for x in range(100):
        ...     fr.run_one_step()
        ...     df.map_depressions()
        ...     flooded = np.where(df.flood_status==3)[0]
        ...     fsc.run_one_step(dt = fsc_dt)
        ...     mg.at_node['topographic__elevation'][0] -= 0.001 #uplift

        Instantiate the E/D component:

        >>> ed = ErosionDeposition(
        ...     mg,
        ...     K=0.00001,
        ...     phi=0.0,
        ...     v_s=0.001,
        ...     m_sp=0.5,
        ...     n_sp = 1.0,
        ...     sp_crit=0,
        ...     erode_flooded_nodes=False)

        Now run the E/D component for 2000 short timesteps:

        >>> for x in range(2000): #E/D component loop
        ...     fr.run_one_step()
        ...     df.map_depressions()
        ...     ed.run_one_step(dt = ed_dt)
        ...     mg.at_node['topographic__elevation'][0] -= 2e-4 * ed_dt

        Now we test to see if topography is right:

        >>> np.around(mg.at_node['topographic__elevation'], decimals=3) # doctest: +NORMALIZE_WHITESPACE
        array([-0.477,  1.036,  2.073,  3.011,  4.061,  1.082, -0.08 , -0.065,
           -0.054,  5.01 ,  2.04 , -0.065, -0.065, -0.053,  6.026,  3.059,
           -0.054, -0.053, -0.035,  7.053,  4.059,  5.041,  6.07 ,  7.004,
            8.01 ])
        """
        if grid.at_node["flow__receiver_node"].size != grid.size("node"):
            msg = (
                "A route-to-multiple flow director has been "
                "run on this grid. The landlab development team has not "
                "verified that ErosionDeposition is compatible with "
                "route-to-multiple methods. Please open a GitHub Issue "
                "to start this process."
            )
            raise NotImplementedError(msg)

        super(ErosionDeposition, self).__init__(
            grid,
            m_sp=m_sp,
            n_sp=n_sp,
            phi=phi,
            F_f=F_f,
            v_s=v_s,
            dt_min=dt_min,
            discharge_field=discharge_field,
            erode_flooded_nodes=erode_flooded_nodes,
        )

        # E/D specific inits.

        # K's and critical values can be floats, grid fields, or arrays
        self._K = return_array_at_node(grid, K)
        self._sp_crit = return_array_at_node(grid, sp_crit)

        # Handle option for solver
        if solver == "basic":
            self.run_one_step = self.run_one_step_basic
        elif solver == "adaptive":
            self.run_one_step = self.run_with_adaptive_time_step_solver
            self._time_to_flat = np.zeros(grid.number_of_nodes)
        else:
            raise ValueError(
                "Parameter 'solver' must be one of: " + "'basic', 'adaptive'"
            )

    def _calc_erosion_rates(self):
        """Calculate erosion rates."""
        omega = self._K * self._Q_to_the_m * np.power(self._slope, self._n_sp)
        omega_over_sp_crit = np.divide(
            omega, self._sp_crit, out=np.zeros_like(omega), where=self._sp_crit != 0
        )

        self._erosion_term = omega - self._sp_crit * (1.0 - np.exp(-omega_over_sp_crit))

    def run_one_step_basic(self, dt=1.0):
        """Calculate change in rock and alluvium thickness for a time period
        'dt'.

        Parameters
        ----------
        dt : float
            Model timestep [T]
        """

        self._calc_hydrology()
        self._calc_erosion_rates()

        if not self._erode_flooded_nodes:
            flood_status = self._grid.at_node["flood_status_code"]
            flooded_nodes = np.nonzero(flood_status == _FLOODED)[0]
        else:
            flooded_nodes = []

        self._erosion_term[flooded_nodes] = 0.0
        self._qs_in[:] = 0.0

        # iterate top to bottom through the stack, calculate qs
        # cythonized version of calculating qs_in
        calculate_qs_in(
            np.flipud(self._stack),
            self._flow_receivers,
            self._cell_area_at_node,
            self._q,
            self._qs,
            self._qs_in,
            self._erosion_term,
            self._v_s,
            self._F_f,
        )

        self._depo_rate[:] = 0.0
        self._depo_rate[self._q > 0] = self._qs[self._q > 0] * (
            self._v_s / self._q[self._q > 0]
        )

        # topo elev is old elev + deposition - erosion
        cores = self._grid.core_nodes
        self._topographic__elevation[cores] += (
            (self._depo_rate[cores] - self._erosion_term[cores]) / (1 - self._phi)
        ) * dt

    def run_with_adaptive_time_step_solver(self, dt=1.0):
        """CHILD-like solver that adjusts time steps to prevent slope
        flattening."""

        # Initialize remaining_time, which records how much of the global time
        # step we have yet to use up.
        remaining_time = dt

        z = self._grid.at_node["topographic__elevation"]
        r = self._flow_receivers
        dzdt = np.zeros(len(z))
        cores = self._grid.core_nodes

        first_iteration = True

        if not self._erode_flooded_nodes:
            flood_status = self._grid.at_node["flood_status_code"]
            flooded_nodes = np.nonzero(flood_status == _FLOODED)[0]
        else:
            flooded_nodes = []

        # Outer WHILE loop: keep going until time is used up
        while remaining_time > 0.0:

            # Update all the flow-link slopes.
            #
            # For the first iteration, we assume this has already been done
            # outside the component (e.g., by flow router), but we need to do
            # it ourselves on subsequent iterations.
            if not first_iteration:
                # update the link slopes
                self._update_flow_link_slopes()
                # update where nodes are flooded. This shouuldn't happen bc
                # of the dynamic timestepper, but just incase, we update here.
                new_flooded_nodes = np.where(self._slope < 0)[0]
                flooded_nodes = np.asarray(
                    np.unique(np.concatenate((flooded_nodes, new_flooded_nodes))),
                    dtype=np.int64,
                )
            else:
                first_iteration = False

            # Calculate rates of entrainment
            self._calc_hydrology()
            self._calc_erosion_rates()
            self._erosion_term[flooded_nodes] = 0.0
            self._qs_in[:] = 0.0

            # Sweep through nodes from upstream to downstream, calculating Qs.
            calculate_qs_in(
                np.flipud(self._stack),
                self._flow_receivers,
                self._cell_area_at_node,
                self._q,
                self._qs,
                self._qs_in,
                self._erosion_term,
                self._v_s,
                self._F_f,
            )

            # Use Qs to calculate deposition rate at each node.
            self._depo_rate[:] = 0.0
            self._depo_rate[self._q > 0] = self._qs[self._q > 0] * (
                self._v_s / self._q[self._q > 0]
            )

            # Rate of change of elevation at core nodes:
            dzdt[cores] = (self._depo_rate[cores] - self._erosion_term[cores]) / (
                1 - self._phi
            )

            # Difference in elevation between each upstream-downstream pair
            zdif = z - z[r]

            # Rate of change of the *difference* in elevation between each
            # upstream-downstream pair.
            rocdif = dzdt - dzdt[r]

            # (Re)-initialize the array that will contain "time to (almost)
            # flat" for each node (relative to its downstream neighbor).
            self._time_to_flat[:] = remaining_time

            # Find locations where the upstream and downstream node elevations
            # are converging (e.g., the upstream one is eroding faster than its
            # downstream neighbor)
            converging = np.where(rocdif < 0.0)[0]

            # Find the time to (almost) flat by dividing difference by rate of
            # change of difference, and then multiplying by a "safety factor"
            self._time_to_flat[converging] = -(
                TIME_STEP_FACTOR * zdif[converging] / rocdif[converging]
            )

            # Mask out pairs where the source at the same or lower elevation
            # as its downstream neighbor (e.g., because it's a pit or a lake).
            # Here, masking out means simply assigning the remaining time in
            # the global time step.
            self._time_to_flat[np.where(zdif <= 0.0)[0]] = remaining_time
            self._time_to_flat[flooded_nodes] = remaining_time

            # From this, find the maximum stable time step. If it is smaller
            # than our tolerance, report and quit.
            dt_max = max(np.amin(self._time_to_flat), self._dt_min)

            # Finally, apply dzdt to all nodes for a (sub)step of duration
            # dt_max
            z[cores] += dzdt[cores] * dt_max

            # Update remaining time and continue the loop
            remaining_time -= dt_max

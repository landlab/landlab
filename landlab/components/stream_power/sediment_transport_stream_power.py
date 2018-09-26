import numpy as np
from landlab.components.erosion_deposition.generalized_erosion_deposition import (_GeneralizedErosionDeposition,
                                                            DEFAULT_MINIMUM_TIME_STEP)
from landlab.utils.return_array import return_array_at_node
from .cfuncs import calculate_qs_in

ROOT2 = np.sqrt(2.0)    # syntactic sugar for precalculated square root of 2
TIME_STEP_FACTOR = 0.5  # factor used in simple subdivision solver
DEFAULT_MINIMUM_TIME_STEP = 0.001  # default minimum time step duration


class TransportLimitedEroder(_GeneralizedErosionDeposition):
    """
    Implementation of a river erosion component which treats fluvial sediment
    flux as a function of stream power, and then uses mass balance (i.e.,
    the divergence of sediment flux) to dictate the response of the river bed.

    In other words,

    Qs = K * A**m * S**n - optional_threshold
    dz/dt = d/dx(Qs)

    If provided, the threshold is implemented smoothly. See the
    StreamPowerSmoothThresholdEroder for more details.
    """
    _name = 'TransportLimitedEroder'

    _input_var_names = (
        'flow__receiver_node',
        'flow__upstream_node_order',
        'topographic__steepest_slope',
        'drainage_area',
    )

    _output_var_names = (
        'topographic__elevation'
    )

    _var_units = {
        'flow__receiver_node': '-',
        'flow__upstream_node_order': '-',
        'topographic__steepest_slope': '-',
        'drainage_area': 'm**2',
        'topographic__elevation': 'm',
    }

    _var_mapping = {
        'flow__receiver_node': 'node',
        'flow__upstream_node_order': 'node',
        'topographic__steepest_slope': 'node',
        'drainage_area': 'node',
        'topographic__elevation': 'node',
    }

    _var_doc = {
        'flow__receiver_node':
            'Node array of receivers (node that receives flow from current '
            'node)',
        'flow__upstream_node_order':
            'Node array containing downstream-to-upstream ordered list of '
            'node IDs',
        'topographic__steepest_slope':
            'Topographic slope at each node',
        'drainage_area':
            "Upstream accumulated surface area contributing to the node's "
            "discharge",
        'topographic__elevation':
            'Land surface topographic elevation',
    }

    def __init__(self, grid, K=None, phi=None, v_s=None,
                 m_sp=None, n_sp=None, sp_crit=0.0, F_f=0.0,
                 discharge_field=None,
                 solver='adaptive',
                 dt_min=DEFAULT_MINIMUM_TIME_STEP,
                 **kwds):
        """Initialize the TransportLimitedEroder component.

        Parameters
        ----------
        grid : ModelGrid
            Landlab ModelGrid object
        K : float, field name, or array
            Prefactor for sediment transport equation (units vary).
        phi : float
            Sediment porosity [-].
        m_sp : float
            Drainage area exponent for sediment transport equation [-]
        n_sp : float
            Slope exponent for sediment transport equation [-]
        sp_crit : float, field name, or array
            Critical stream power to erode substrate [E/(TL^2)]
        F_f : float
            Fraction of eroded material that turns into "fines" that do not
            contribute to (coarse) sediment load. Defaults to zero.
        discharge_field : float, field name, or array, or None
            If None, the component will use its standard equation for sediment
            transport, Qc = K * A**m * S**n. If it has a value, instead the
            component will implement Qc = K * Q**m * S**n, where Q is
            discharge_field, a water discharge [L^2/T].
        solver : string
            Solver to use. Options at present include:
                (1) 'basic': explicit forward-time extrapolation.
                    Simple but will become unstable if time step is too large.
                (2) 'adaptive' (default): adaptive time-step solver that
                    estimates a stable step size based on the shortest time to
                    "flattening" among all upstream-downstream node pairs.
        """
        if (grid.at_node['flow__receiver_node'].size != grid.size('node')):
            msg = ('A route-to-multiple flow director has been '
                   'run on this grid. The landlab development team has not '
                   'verified that ErosionDeposition is compatible with '
                   'route-to-multiple methods. Please open a GitHub Issue '
                   'to start this process.')
            raise NotImplementedError(msg)

        if discharge_field is not None:
            super(ErosionDeposition, self).__init__(
                grid, m_sp=m_sp, n_sp=n_sp, phi=phi, F_f=F_f, v_s=1.,
                dt_min=dt_min, discharge_field=discharge_field)
        else:
            super(ErosionDeposition, self).__init__(
                grid, m_sp=m_sp, n_sp=n_sp, phi=phi, F_f=F_f, v_s=1.,
                dt_min=dt_min, discharge_field='drainage_area')

        self._grid = grid  # store grid

        # K's and critical values can be floats, grid fields, or arrays
        self.K = return_array_at_node(grid, K)
        # special cases if sp_crit is 0 or <0...
        self.sp_crit = return_array_at_node(grid, sp_crit)
        if np.any(self.sp_crit < 0.):
            raise ValueError('sp_crit must be >= 0. everywhere')
        if np.allclose(sp_crit, 0.):
            self._calc_trp = self._calc_transport_rates_wo_thresh()
        else:
            self._calc_trp = self._calc_transport_rates_with_thresh()

        # Handle option for solver
        if solver == 'basic':
            self.run_one_step = self.run_one_step_basic
        elif solver == 'adaptive':
            self.run_one_step = self.run_with_adaptive_time_step_solver
            self.time_to_flat = np.zeros(grid.number_of_nodes)
        else:
            raise ValueError("Parameter 'solver' must be one of: "
                             + "'basic', 'adaptive'")

    def _calc_transport_rates_with_thresh(self):
        """Calculate Qs with a threshold"""
        omega = self.K * self.Q_to_the_m * np.power(self.slope, self.n_sp)
        omega_over_sp_crit = np.divide(
            omega, self.sp_crit, out=np.zeros_like(omega),
            where=np.logical_not(np.isclose(self.sp_crit, 0.)))

        self.qs = omega - self.sp_crit * (1.0 - np.exp(-omega_over_sp_crit))

    def _calc_transport_rates_wo_thresh(self):
        """Calculate Qs with a threshold"""
        self.qs = self.K * self.Q_to_the_m * np.power(self.slope, self.n_sp)

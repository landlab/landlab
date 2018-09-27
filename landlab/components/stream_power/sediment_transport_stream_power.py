import numpy as np
from landlab.components.erosion_deposition.generalized_erosion_deposition import (_GeneralizedErosionDeposition,
                                                            DEFAULT_MINIMUM_TIME_STEP)
from landlab.utils.return_array import return_array_at_node

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

    The component permits supply of a porosity (phi) and wash fraction (F_f).
    Note that if these are not 0., then you are implicitly assuming that
    any sediment pickup is occurring out of intact rock, that has not
    previously been processed by a river. This assumption may often not be
    appropriate.
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

    def __init__(self, grid, K=None,
                 m_sp=0.5, n_sp=1., sp_crit=0.0,  phi=0., F_f=0.0,
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
            super(TransportLimitedEroder, self).__init__(
                grid, m_sp=m_sp, n_sp=n_sp, phi=phi, F_f=F_f, v_s=1.,
                dt_min=dt_min, discharge_field=discharge_field)
        else:
            super(TransportLimitedEroder, self).__init__(
                grid, m_sp=m_sp, n_sp=n_sp, phi=phi, F_f=F_f, v_s=1.,
                dt_min=dt_min, discharge_field='drainage_area')

        self._grid = grid  # store grid
        sed_fract_to_remain = (1. - F_f) * (1. - phi)
        self._one_by_erosion_loss = 1. / sed_fract_to_remain
        # tests in _GeneralizedErosionDeposition stop this ever exploding
        if np.isclose(self._one_by_erosion_loss, 1.):
            self._calc_sed_div = _calc_sed_flux_divergence
        else:
            self._calc_sed_div = _calc_sed_flux_divergence_lossy

        # the net loss of rock during erosion

        # K's and critical values can be floats, grid fields, or arrays
        self.K = return_array_at_node(grid, K)
        self.Qs_in = grid.zeros('node', dtype=float)
        # special cases if sp_crit is 0 or <0...
        self.sp_crit = return_array_at_node(grid, sp_crit)
        if np.any(self.sp_crit < 0.):
            raise ValueError('sp_crit must be >= 0. everywhere')
        if np.allclose(sp_crit, 0.):
            self._calc_trp = self._calc_transport_rates_wo_thresh
        else:
            self._calc_trp = self._calc_transport_rates_with_thresh

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

    def _calc_discharges(self, flooded_nodes):
        """
        Calculate sediment discharge out of a node in channel.
        """
        self._calc_trp()  # this sets fluxes, so need to weight by A
        self.Qs_out = self.qs * self.cell_area_at_node
        if flooded_nodes is not None:
            self.Qs_out[flooded_nodes] = 0.

    def run_with_adaptive_time_step_solver(self, dt=1.0, flooded_nodes=[],
                                           **kwds):
        """CHILD-like solver that adjusts time steps to prevent slope
        flattening."""

        # Initialize remaining_time, which records how much of the global time
        # step we have yet to use up.
        remaining_time = dt

        # z = self._grid.at_node['topographic__elevation']
        r = self.flow_receivers
        dzdt = self.grid.zeros('node', dtype=float)
        cores = self.grid.core_nodes

        first_iteration = True

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
                new_flooded_nodes = np.where(self.slope < 0)[0]
                if len(new_flooded_nodes != 0):
                    flooded_nodes = np.asarray(np.unique(
                        np.concatenate((flooded_nodes, new_flooded_nodes))),
                        dtype=np.int64)
            else:
                first_iteration = False

            self.Qs_in[:] = 0.0

            self._calc_discharges(flooded_nodes=flooded_nodes)
            self._calc_sed_div(np.flipud(self.stack), r, self.Qs_out,
                               self.Qs_in, self._one_by_erosion_loss)

            # Qs_in now gives dz/dt:
            dzdt[cores] = (
                self.Qs_in[cores] / self.grid.cell_area_at_node[cores])

            # Difference in elevation between each upstream-downstream pair
            zdif = z - z[r]
            # Rate of change of the *difference* in elevation between each
            # upstream-downstream pair.
            rocdif = dzdt - dzdt[r]

            # (Re)-initialize the array that will contain "time to (almost)
            # flat" for each node (relative to its downstream neighbor).
            self.time_to_flat[:] = remaining_time

            # Find locations where the upstream and downstream node elevations
            # are converging (e.g., the upstream one is eroding faster than its
            # downstream neighbor)
            converging = rocdif < 0.0

            # Find the time to (almost) flat by dividing difference by rate of
            # change of difference, and then multiplying by a "safety factor"
            self.time_to_flat[converging] = -1. * (TIME_STEP_FACTOR
                                                   * zdif[converging]
                                                   / rocdif[converging])

            # Mask out pairs where the source at the same or lower elevation
            # as its downstream neighbor (e.g., because it's a pit or a lake).
            # Here, masking out means simply assigning the remaining time in
            # the global time step.
            self.time_to_flat[zdif <= 0.0] = remaining_time
            self.time_to_flat[flooded_nodes] = remaining_time

            # From this, find the maximum stable time step. If it is smaller
            # than our tolerance, reduce to tolerance
            dt_max = np.amin(self.time_to_flat)
            if dt_max < self.dt_min:
                dt_max = self.dt_min

            # TODO: this is CMS's method, and it's not terribly efficient.
            # Could enhance to remove need to search whole grid.

            # Finally, apply dzdt to all nodes for a (sub)step of duration
            # dt_max
            z[cores] += dzdt[cores] * dt_max

            # Update remaining time and continue the loop
            remaining_time -= dt_max


# These following two methods are designed to be cythonised, but written
# like this to enable clean testing.


def _calc_sed_flux_divergence_lossy(stack_up_to_down,
                                    flow_receivers,
                                    Qs,
                                    Qs_in,
                                    one_by_erosion_loss):
    """
    Calculate divergence of sediment discharge at a node in channel.
    Qs_in should begin as either an array of zeros, or potentially an
    existing, non-fluvial sediment discharge into the node.
    This version accounts for loss of material during erosion.
    """
    # cdef unsigned int node_id
    # cdef unsigned int next_node
    # cdef double Qs_here
    # We iterate this calc to correctly add the discharges dstr
    for node_id in stack_up_to_down:
        next_node = flow_receivers[node_id]
        Qs_here = Qs[node_id]
        Qs_in[next_node] += Qs_here
        Qs_in[node_id] -= Qs_here   # note these cancel out in a self-drainer
        if Qs_in[node_id] < 0.:
            # erosion is occurring. Need to consider F_f & phi.
            Qs_in[node_id] *= one_by_erosion_loss


def _calc_sed_flux_divergence(stack_up_to_down,
                              flow_receivers,
                              Qs,
                              Qs_in,
                              dummy):
    """
    Calculate divergence of sediment discharge at a node in channel.
    Qs_in should begin as either an array of zeros, or potentially an
    existing, non-fluvial sediment discharge into the node.
    """
    # cdef unsigned int node_id
    # cdef unsigned int next_node
    # cdef double Qs_here
    # We iterate this calc to correctly add the discharges dstr
    for node_id in stack_up_to_down:
        next_node = flow_receivers[node_id]
        Qs_here = Qs[node_id]
        Qs_in[next_node] += Qs_here
        Qs_in[node_id] -= Qs_here  # note these cancel out in a self-drainer

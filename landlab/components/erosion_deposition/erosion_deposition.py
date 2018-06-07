import numpy as np
from landlab import Component
from landlab import RasterModelGrid
from .cfuncs import calculate_qs_in

ROOT2 = np.sqrt(2.0)    # syntactic sugar for precalculated square root of 2
TIME_STEP_FACTOR = 0.5  # factor used in simple subdivision solver
DEFAULT_MINIMUM_TIME_STEP = 0.001  # default minimum time step duration

class ErosionDeposition(Component):
    """
    Erosion-Deposition model in the style of Davy and Lague (2009)

    Component written by C. Shobe, begun July 2016.
    """

    _name= 'ErosionDeposition'

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
                 m_sp=None, n_sp=None, sp_crit=None, F_f=0.0,
                 method=None, discharge_method=None,
                 area_field=None, discharge_field=None, solver='basic',
                 dt_min=DEFAULT_MINIMUM_TIME_STEP,
                 **kwds):
        """Initialize the ErosionDeposition model.

        Parameters
        ----------
        grid : ModelGrid
            Landlab ModelGrid object
        K : float
            Erodibility constant for substrate (units vary).
        phi : float
            Sediment porosity [-].
        v_s : float
            Effective settling velocity for chosen grain size metric [L/T].
        m_sp : float
            Drainage area exponent (units vary)
        n_sp : float
            Slope exponent (units vary)
        sp_crit : float
            Critical stream power to erode substrate [E/(TL^2)]
        F_f : float
            Fraction of eroded material that turns into "fines" that do not
            contribute to (coarse) sediment load. Defaults to zero.
        method : string
            Either "simple_stream_power", "threshold_stream_power", or
            "stochastic_hydrology". Method for calculating sediment
            and bedrock entrainment/erosion.
        discharge_method : string
            Either "area_field" or "discharge_field". If using stochastic
            hydrology, determines whether component is supplied with
            drainage area or discharge.
        area_field : string or array
            Used if discharge_method = 'area_field'. Either field name or
            array of length(number_of_nodes) containing drainage areas [L^2].
        discharge_field : string or array
            Used if discharge_method = 'discharge_field'.Either field name or
            array of length(number_of_nodes) containing drainage areas [L^2/T].
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
        >>> from landlab.components.flow_routing import FlowRouter
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
        >>> mg = RasterModelGrid((nr, nc), 10.0)
        >>> _ = mg.add_zeros('node', 'topographic__elevation')
        >>> mg['node']['topographic__elevation'] += mg.node_y/10 + \
                mg.node_x/10 + np.random.rand(len(mg.node_y)) / 10
        >>> mg.set_closed_boundaries_at_grid_edges(bottom_is_closed=True,\
                                                       left_is_closed=True,\
                                                       right_is_closed=True,\
                                                       top_is_closed=True)
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

        >>> fsc = FastscapeEroder(mg, K_sp=.001, m_sp=.5, n_sp=1)
        >>> fr = FlowRouter(mg) #instantiate
        >>> df = DepressionFinderAndRouter(mg)

        Burn in an initial drainage network using the Fastscape eroder:

        >>> for x in range(100):
        ...     fr.run_one_step()
        ...     df.map_depressions()
        ...     flooded = np.where(df.flood_status==3)[0]
        ...     fsc.run_one_step(dt = fsc_dt, flooded_nodes=flooded)
        ...     mg.at_node['topographic__elevation'][0] -= 0.001 #uplift

        Instantiate the E/D component:

        >>> ed = ErosionDeposition(mg, K=0.00001, phi=0.0, v_s=0.001,\
                                m_sp=0.5, n_sp = 1.0, sp_crit=0,\
                                method='simple_stream_power',\
                                discharge_method=None, area_field=None,\
                                discharge_field=None)

        Now run the E/D component for 2000 short timesteps:

        >>> for x in range(2000): #E/D component loop
        ...     fr.run_one_step()
        ...     df.map_depressions()
        ...     flooded = np.where(df.flood_status==3)[0]
        ...     ed.run_one_step(dt = ed_dt, flooded_nodes=flooded)
        ...     mg.at_node['topographic__elevation'][0] -= 2e-4 * ed_dt

        Now we test to see if topography is right:

        >>> np.around(mg.at_node['topographic__elevation'], decimals=3) # doctest: +NORMALIZE_WHITESPACE
        array([-0.477,  1.036,  2.073,  3.011,  4.061,  1.082, -0.08 , -0.065,
           -0.054,  5.01 ,  2.04 , -0.065, -0.065, -0.053,  6.026,  3.059,
           -0.054, -0.053, -0.035,  7.053,  4.059,  5.041,  6.07 ,  7.004,
            8.01 ])
        """
#        array([-0.47709402,  1.03606698,  2.0727653 ,  3.01126678,  4.06077707,
#            1.08157495, -0.0799798 , -0.06459322, -0.05380581,  5.00969486,
#            2.04008677, -0.06457996, -0.06457219, -0.05266169,  6.02641123,
#            3.05874171, -0.05350698, -0.05265586, -0.03498794,  7.05334077,
#            4.05922478,  5.0409473 ,  6.07035008,  7.0038935 ,  8.01034357])
        # assign class variables to grid fields; create necessary fields
        self.flow_receivers = grid.at_node['flow__receiver_node']
        self.stack = grid.at_node['flow__upstream_node_order']
        self.elev = grid.at_node['topographic__elevation']
        self.slope = grid.at_node['topographic__steepest_slope']
        self.link_to_reciever = grid.at_node['flow__link_to_receiver_node']
        self.cell_area_at_node = grid.cell_area_at_node

        if isinstance(grid, RasterModelGrid):
            self.link_lengths = grid.length_of_d8
        else:
            self.link_lengths = grid.length_of_link

        try:
            self.qs = grid.at_node['sediment__flux']
        except KeyError:
            self.qs = grid.add_zeros(
                'sediment__flux', at='node', dtype=float)
        try:
            self.q = grid.at_node['surface_water__discharge']
        except KeyError:
            self.q = grid.add_zeros(
                'surface_water__discharge', at='node', dtype=float)

        self._grid = grid #store grid

        # Create arrays for sediment influx at each node, discharge to the
        # power "m", and deposition rate
        self.qs_in = np.zeros(grid.number_of_nodes)
        self.Q_to_the_m = np.zeros(grid.number_of_nodes)
        self.S_to_the_n = np.zeros(grid.number_of_nodes)
        self.depo_rate = np.zeros(self.grid.number_of_nodes)

        # store other constants
        self.m_sp = float(m_sp)
        self.n_sp = float(n_sp)
        self.phi = float(phi)
        self.v_s = float(v_s)
        self.dt_min = dt_min
        self.frac_coarse = 1.0 - F_f

        # K's and critical values can be floats, grid fields, or arrays
        if type(K) is str:
            self.K = self._grid.at_node[K]
        elif type(K) in (float, int):  # a float
            self.K = float(K)
        elif len(K) == self.grid.number_of_nodes:
            self.K = np.array(K)
        else:
            raise TypeError('Supplied type of K ' +
                            'was not recognised, or array was ' +
                            'not nnodes long!')

        if sp_crit is not None:
            if type(sp_crit) is str:
                self.sp_crit = self._grid.at_node[sp_crit]
            elif type(sp_crit) in (float, int):  # a float
                self.sp_crit = float(sp_crit)
            elif len(sp_crit) == self.grid.number_of_nodes:
                self.sp_crit = np.array(sp_crit)
            else:
                raise TypeError('Supplied type of sp_crit ' +
                                'was not recognised, or array was ' +
                                'not nnodes long!')

        #go through erosion methods to ensure correct hydrology
        self.method = str(method)
        if discharge_method is not None:
            self.discharge_method = str(discharge_method)
        else:
            self.discharge_method = None
        if area_field is not None:
            self.area_field = str(area_field)
        else:
            self.area_field = None
        if discharge_field is not None:
            self.discharge_field = str(discharge_field)
        else:
            self.discharge_field = None

        if self.method == 'simple_stream_power':
            self.calc_ero_rate = self.simple_stream_power
        elif self.method == 'threshold_stream_power':
            self.calc_ero_rate = self.threshold_stream_power
        elif self.method == 'stochastic_hydrology':
            self.calc_ero_rate = self.stochastic_hydrology
        else:
            print('METHOD:')
            print(self.method)
            raise ValueError('Specify erosion method (simple stream power,\
                            threshold stream power, or stochastic hydrology)!')

        # Handle option for solver
        if solver == 'basic':
            self.run_one_step = self.run_one_step_basic
        elif solver == 'adaptive':
            self.run_one_step = self.run_with_adaptive_time_step_solver
            self.time_to_flat = np.zeros(grid.number_of_nodes)
        else:
            raise ValueError("Parameter 'solver' must be one of: "
                             + "'basic', 'adaptive'")

    #three choices for erosion methods:
    def simple_stream_power(self):
        """Use non-threshold stream power.

        simple_stream_power uses no entrainment or erosion thresholds,
        and uses either q=A^m or q=Q^m depending on discharge method. If
        discharge method is None, default is q=A^m.
        """
        #self.Q_to_the_m = np.zeros(len(self.grid.at_node['drainage_area']))
        if self.method == 'simple_stream_power' and self.discharge_method == None:
            self.Q_to_the_m[:] = np.power(self.grid.at_node['drainage_area'], self.m_sp)
        elif self.method == 'simple_stream_power' and self.discharge_method is not None:
            if self.discharge_method == 'drainage_area':
                if self.area_field is not None:
                    if type(self.area_field) is str:
                        self.drainage_area = self._grid.at_node[self.area_field]
                    elif len(self.area_field) == self.grid.number_of_nodes:
                        self.drainage_area = np.array(self.area_field)
                    else:
                        raise TypeError('Supplied type of area_field ' +
                                'was not recognised, or array was ' +
                                'not nnodes long!')
                self.Q_to_the_m[:] = np.power(self.drainage_area, self.m_sp)
            elif self.discharge_method == 'discharge_field':
                if self.discharge_field is not None:
                    if type(self.discharge_field) is str:
                        self.q[:] = self._grid.at_node[self.discharge_field]
                        self.Q_to_the_m[:] = np.power(self.q, self.m_sp)
                    elif len(self.discharge_field) == self.grid.number_of_nodes:
                        self.q[:] = np.array(self.discharge_field)
                        self.Q_to_the_m[:] = np.power(self.q, self.m_sp)
                    else:
                        raise TypeError('Supplied type of discharge_field ' +
                                'was not recognised, or array was ' +
                                'not nnodes long!')
        self.S_to_the_n[:] = 0
        self.S_to_the_n[self.slope > 0] = np.power(self.slope[self.slope > 0] , self.n_sp)
        self.erosion_term = self.K * self.Q_to_the_m * self.S_to_the_n

        self.qs_in[:] = 0.0

    def threshold_stream_power(self):
        """Use stream power with entrainment/erosion thresholds.

        threshold_stream_power works the same way as simple SP but includes
        user-defined thresholds for sediment entrainment and bedrock erosion.
        """
        #self.Q_to_the_m = np.zeros(len(self.grid.at_node['drainage_area']))
        if self.method == 'threshold_stream_power' and self.discharge_method == None:
            self.Q_to_the_m[:] = np.power(self.grid.at_node['drainage_area'], self.m_sp)
        elif self.method == 'threshold_stream_power' and self.discharge_method is not None:
            if self.discharge_method == 'drainage_area':
                if self.area_field is not None:
                    if type(self.area_field) is str:
                        self.drainage_area = self._grid.at_node[self.area_field]
                    elif len(self.area_field) == self.grid.number_of_nodes:
                        self.drainage_area = np.array(self.area_field)
                    else:
                        raise TypeError('Supplied type of area_field ' +
                                'was not recognised, or array was ' +
                                'not nnodes long!')
                self.Q_to_the_m[:] = np.power(self.drainage_area, self.m_sp)
            elif self.discharge_method == 'discharge_field':
                if self.discharge_field is not None:
                    if type(self.discharge_field) is str:
                        self.q[:] = self._grid.at_node[self.discharge_field]
                        self.Q_to_the_m[:] = np.power(self.q, self.m_sp)
                    elif len(self.discharge_field) == self.grid.number_of_nodes:
                        self.q[:] = np.array(self.discharge_field)
                        self.Q_to_the_m[:] = np.power(self.q, self.m_sp)
                    else:
                        raise TypeError('Supplied type of discharge_field ' +
                                'was not recognised, or array was ' +
                                'not nnodes long!')

        self.S_to_the_n[:] = 0
        self.S_to_the_n[self.slope > 0] = np.power(self.slope[self.slope > 0] , self.n_sp)
        self.erosion_term = self.K * self.Q_to_the_m * self.S_to_the_n

        omega = self.K * self.Q_to_the_m * self.S_to_the_n

        self.erosion_term = omega - self.sp_crit * \
            (1 - np.exp(-omega / self.sp_crit))
        self.qs_in[:] = 0.0

    def stochastic_hydrology(self):
        """Allows custom area and discharge fields, no default behavior.

        stochastic_hydrology forces the user to supply either an array or
        field name for either drainage area or discharge, and will not
        default to q=A^m.
        """
        if self.method == 'stochastic_hydrology' and self.discharge_method == None:
            raise TypeError('Supply a discharge method to use stoc. hydro!')
        elif self.discharge_method is not None:
            if self.discharge_method == 'drainage_area':
                if self.area_field is not None:
                    if type(self.area_field) is str:
                        self.drainage_area = self._grid.at_node[self.area_field]
                    elif len(self.area_field) == self.grid.number_of_nodes:
                        self.drainage_area = np.array(self.area_field)
                    else:
                        raise TypeError('Supplied type of area_field ' +
                                'was not recognised, or array was ' +
                                'not nnodes long!')
                self.Q_to_the_m[:] = np.power(self.grid.at_node['drainage_area'], self.m_sp)
            elif self.discharge_method == 'discharge_field':
                if self.discharge_field is not None:
                    if type(self.discharge_field) is str:
                        self.q[:] = self._grid.at_node[self.discharge_field]
                        self.Q_to_the_m[:] = np.power(self.q, self.m_sp)
                    elif len(self.discharge_field) == self.grid.number_of_nodes:
                        self.q[:] = np.array(self.discharge_field)
                        self.Q_to_the_m[:] = np.power(self.q, self.m_sp)
                    else:
                        raise TypeError('Supplied type of discharge_field ' +
                                'was not recognised, or array was ' +
                                'not nnodes long!')
            else:
                raise ValueError('Specify discharge method for stoch hydro!')
        self.S_to_the_n[:] = 0
        self.S_to_the_n[self.slope > 0] = np.power(self.slope[self.slope > 0] , self.n_sp)
        self.erosion_term = self.K * self.Q_to_the_m * self.S_to_the_n

        self.qs_in[:] = 0.0

    def _update_flow_link_slopes(self):
        """Updates gradient between each core node and its receiver.

        Used to update slope values between sub-time-steps, when we do not
        re-run flow routing.

        Examples
        --------
        >>> from landlab import RasterModelGrid
        >>> from landlab.components import FlowAccumulator
        >>> rg = RasterModelGrid((3, 4))
        >>> z = rg.add_zeros('node', 'topographic__elevation')
        >>> z[:] = rg.x_of_node + rg.y_of_node
        >>> fa = FlowAccumulator(rg, flow_director='FlowDirectorD8')
        >>> fa.run_one_step()
        >>> rg.at_node['topographic__steepest_slope'][5:7]
        array([ 1.41421356,  1.41421356])
        >>> sp = ErosionDeposition(rg, K=0.00001, phi=0.1, v_s=0.001,\
                                   m_sp=0.5, n_sp = 1.0, sp_crit_sed=0,\
                                   sp_crit_br=0, method='simple_stream_power',\
                                   discharge_method=None, area_field=None,\
                                   discharge_field=None)
        >>> z *= 0.1
        >>> sp._update_flow_link_slopes()
        >>> rg.at_node['topographic__steepest_slope'][5:7]
        array([ 0.14142136,  0.14142136])
        """
        z = self._grid.at_node['topographic__elevation']
        r = self._grid.at_node['flow__receiver_node']
        slp = self._grid.at_node['topographic__steepest_slope']
        slp[:] = (z - z[r]) / self.link_lengths[self.link_to_reciever]

    def run_one_step_basic(self, dt=1.0, flooded_nodes=[], **kwds):
        """Calculate change in rock and alluvium thickness for
           a time period 'dt'.

        Parameters
        ----------
        dt : float
            Model timestep [T]
        flooded_nodes : array
            Indices of flooded nodes, passed from flow router
        """

        self.calc_ero_rate()
        self.erosion_term[flooded_nodes] = 0.0
        self.qs_in[:] = 0.0

        #iterate top to bottom through the stack, calculate qs
        # cythonized version of calculating qs_in
        calculate_qs_in(np.flipud(self.stack),
                        self.flow_receivers,
                        self.cell_area_at_node,
                        self.q,
                        self.qs,
                        self.qs_in,
                        self.erosion_term,
                        self.v_s,
                        self.frac_coarse,
                        self.phi)

        self.depo_rate[:] = 0.0
        self.depo_rate[self.q > 0] = (self.qs[self.q > 0] * \
                                         (self.v_s / self.q[self.q > 0]))

        #topo elev is old elev + deposition - erosion
        cores = self.grid.core_nodes
        self.elev[cores] += ((self.depo_rate[cores]
                              - self.erosion_term[cores]) * dt)

    def run_with_adaptive_time_step_solver(self, dt=1.0, flooded_nodes=[],
                                           **kwds):
        """CHILD-like solver that adjusts time steps to prevent slope
        flattening."""

        # Initialize remaining_time, which records how much of the global time
        # step we have yet to use up.
        remaining_time = dt

        z = self._grid.at_node['topographic__elevation']
        r = self.flow_receivers
        dzdt = np.zeros(len(z))
        cores = self._grid.core_nodes

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
                # update where nodes are flooded. This shouuldn't happen because
                # of the dynamic timestepper, but just incase, we update here.
                new_flooded_nodes = np.where(self.slope<0)[0]
                flooded_nodes = np.asarray(np.unique(np.concatenate((flooded_nodes,
                                                          new_flooded_nodes))), dtype=np.int64)
            else:
                first_iteration = False

            # Calculate rates of entrainment
            self.calc_ero_rate()
            self.erosion_term[flooded_nodes] = 0.0

            # Sweep through nodes from upstream to downstream, calculating Qs.
            calculate_qs_in(np.flipud(self.stack),
                            self.flow_receivers,
                            self.cell_area_at_node,
                            self.q,
                            self.qs,
                            self.qs_in,
                            self.erosion_term,
                            self.v_s,
                            self.frac_coarse,
                            self.phi)

            # Use Qs to calculate deposition rate at each node.
            self.depo_rate[:] = 0.0
            self.depo_rate[self.q > 0] = (self.qs[self.q > 0]
                                          * (self.v_s / self.q[self.q > 0]))

            # Rate of change of elevation at core nodes:
            dzdt[cores] = self.depo_rate[cores] - self.erosion_term[cores]

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
            converging = np.where(rocdif < 0.0)[0]

            # Find the time to (almost) flat by dividing difference by rate of
            # change of difference, and then multiplying by a "safety factor"
            self.time_to_flat[converging] = - (TIME_STEP_FACTOR
                                               * zdif[converging]
                                              / rocdif[converging])

            # Mask out pairs where the source at the same or lower elevation
            # as its downstream neighbor (e.g., because it's a pit or a lake).
            # Here, masking out means simply assigning the remaining time in
            # the global time step.
            self.time_to_flat[np.where(zdif <= 0.0)[0]] = remaining_time
            self.time_to_flat[flooded_nodes] = remaining_time

            # From this, find the maximum stable time step. If it is smaller
            # than our tolerance, report and quit.
            dt_max = np.amin(self.time_to_flat)
            if dt_max < self.dt_min:
                dt_max = self.dt_min

            # Finally, apply dzdt to all nodes for a (sub)step of duration
            # dt_max
            z[cores] += dzdt[cores] * dt_max

            # Update remaining time and continue the loop
            remaining_time -= dt_max

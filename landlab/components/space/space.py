from six import string_types

import numpy as np
from landlab import Component
from landlab import RasterModelGrid
from .cfuncs import calculate_qs_in

ROOT2 = np.sqrt(2.0) # syntactic sugar for precalculated square root of 2
TIME_STEP_FACTOR = 0.5 # factor used in simple subdivision solver
DEFAULT_MINIMUM_TIME_STEP = 0.001 # default minimum time step duration


class Space(Component):
    """Stream Power with Alluvium Conservation and Entrainment (SPACE)

    See the publication:
    Shobe, C. M., Tucker, G. E., and Barnhart, K. R.: The SPACE 1.0 model: a
    Landlab component for 2-D calculation of sediment transport, bedrock
    erosion, and landscape evolution, Geosci. Model Dev., 10, 4577-4604,
    https://doi.org/10.5194/gmd-10-4577-2017, 2017.

    Parameters
    ----------
    grid : ModelGrid
        Landlab ModelGrid object
    K_sed : float
        Erodibility constant for sediment (units vary).
    K_br : float
        Erodibility constant for bedrock (units vary).
    F_f : float
        Fraction of permanently suspendable fines in bedrock [-].
    phi : float
        Sediment porosity [-].
    H_star : float
        Sediment thickness required for full entrainment [L].
    v_s : float
        Effective settling velocity for chosen grain size metric [L/T].
    m_sp : float
        Drainage area exponent (units vary)
    n_sp : float
        Slope exponent (units vary)
    sp_crit_sed : float
        Critical stream power to erode sediment [E/(TL^2)]
    sp_crit_br : float
        Critical stream power to erode rock [E/(TL^2)]
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
            (2) 'adaptive': subdivides global time step as needed to
                prevent slopes from reversing and alluvium from going
                negative.

    Examples
    ---------
    >>> import numpy as np
    >>> from landlab import RasterModelGrid
    >>> from landlab.components.flow_routing import FlowRouter
    >>> from landlab.components import DepressionFinderAndRouter
    >>> from landlab.components import Space
    >>> from landlab.components import FastscapeEroder
    >>> np.random.seed(seed = 5000)

    Define grid and initial topography:

    *  5x5 grid with baselevel in the lower left corner
    *  All other boundary nodes closed
    *  Initial topography is plane tilted up to the upper right with
       noise

    >>> mg = RasterModelGrid((5, 5), spacing=10.0)
    >>> _ = mg.add_zeros('topographic__elevation', at='node')
    >>> mg.at_node['topographic__elevation'] += (mg.node_y / 10. +
    ...     mg.node_x / 10. + np.random.rand(len(mg.node_y)) / 10.)
    >>> mg.set_closed_boundaries_at_grid_edges(bottom_is_closed=True,
    ...                                        left_is_closed=True,
    ...                                        right_is_closed=True,
    ...                                        top_is_closed=True)
    >>> mg.set_watershed_boundary_condition_outlet_id(
    ...     0, mg.at_node['topographic__elevation'], -9999.)
    >>> fsc_dt = 100.
    >>> space_dt = 100.

    Instantiate Fastscape eroder, flow router, and depression finder

    >>> fsc = FastscapeEroder(mg, K_sp=.001, m_sp=.5, n_sp=1)
    >>> fr = FlowRouter(mg)
    >>> df = DepressionFinderAndRouter(mg)

    Burn in an initial drainage network using the Fastscape eroder:

    >>> for x in range(100):
    ...     fr.run_one_step()
    ...     df.map_depressions()
    ...     flooded = np.where(df.flood_status == 3)[0]
    ...     fsc.run_one_step(dt=fsc_dt, flooded_nodes=flooded)
    ...     mg.at_node['topographic__elevation'][0] -= 0.001 # Uplift

    Add some soil to the drainage network:

    >>> _ = mg.add_zeros('soil__depth', at='node', dtype=float)
    >>> mg.at_node['soil__depth'] += 0.5
    >>> mg.at_node['topographic__elevation'] += mg.at_node['soil__depth']

    Instantiate the Space component:

    >>> ha = Space(mg, K_sed=0.00001, K_br=0.00000000001,
    ...            F_f=0.5, phi=0.1, H_star=1., v_s=0.001,
    ...            m_sp=0.5, n_sp = 1.0, sp_crit_sed=0,
    ...            sp_crit_br=0, method='simple_stream_power',
    ...            discharge_method=None, area_field=None,
    ...            discharge_field=None)

    Now run the Space component for 2000 short timesteps:

    >>> for x in range(2000): #Space component loop
    ...     fr.run_one_step()
    ...     df.map_depressions()
    ...     flooded = np.where(df.flood_status == 3)[0]
    ...     ha.run_one_step(dt=space_dt, flooded_nodes=flooded)
    ...     mg.at_node['bedrock__elevation'][0] -= 2e-6 * space_dt

    Now we test to see if soil depth and topography are right:

    >>> np.around(mg.at_node['soil__depth'], decimals=3) # doctest: +NORMALIZE_WHITESPACE
    array([ 0.5  ,  0.5  ,  0.5  ,  0.5  ,  0.5  ,  0.5  ,  0.495,  0.493,
            0.492,  0.5  ,  0.5  ,  0.493,  0.493,  0.491,  0.5  ,  0.5  ,
            0.492,  0.491,  0.486,  0.5  ,  0.5  ,  0.5  ,  0.5  ,  0.5  ,  0.5  ])

    >>> np.around(mg.at_node['topographic__elevation'], decimals=3) # doctest: +NORMALIZE_WHITESPACE
    array([ 0.423,  1.536,  2.573,  3.511,  4.561,  1.582,  0.424,  0.429,
            0.438,  5.51 ,  2.54 ,  0.429,  0.429,  0.439,  6.526,  3.559,
            0.438,  0.439,  0.451,  7.553,  4.559,  5.541,  6.57 ,  7.504,
            8.51 ])
    """

    _name= 'Space'

    _input_var_names = (
        'flow__receiver_node',
        'flow__upstream_node_order',
        'topographic__steepest_slope',
        'drainage_area',
        'soil__depth'
    )

    _output_var_names = (
        'topographic__elevation'
        'soil__depth'
    )

    _var_units = {
        'flow__receiver_node': '-',
        'flow__upstream_node_order': '-',
        'topographic__steepest_slope': '-',
        'drainage_area': 'm**2',
        'soil__depth': 'm',
        'topographic__elevation': 'm',
    }

    _var_mapping = {
        'flow__receiver_node': 'node',
        'flow__upstream_node_order': 'node',
        'topographic__steepest_slope': 'node',
        'drainage_area': 'node',
        'soil__depth': 'node',
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
        'soil__depth':
            'Depth of sediment above bedrock',
        'topographic__elevation':
            'Land surface topographic elevation',
    }

    _cite_as = """@Article{gmd-10-4577-2017,
                  AUTHOR = {Shobe, C. M. and Tucker, G. E. and Barnhart, K. R.},
                  TITLE = {The SPACE~1.0 model: a~Landlab component for 2-D calculation of sediment transport, bedrock erosion, and landscape evolution},
                  JOURNAL = {Geoscientific Model Development},
                  VOLUME = {10},
                  YEAR = {2017},
                  NUMBER = {12},
                  PAGES = {4577--4604},
                  URL = {https://www.geosci-model-dev.net/10/4577/2017/},
                  DOI = {10.5194/gmd-10-4577-2017}
                  }"""

    def __init__(self, grid, K_sed=None, K_br=None, F_f=None,
                 phi=None, H_star=None, v_s=None,
                 m_sp=None, n_sp=None, sp_crit_sed=None,
                 sp_crit_br=None, method=None, discharge_method=None,
                 area_field=None, discharge_field=None, solver='basic',
                 dt_min=DEFAULT_MINIMUM_TIME_STEP,
                 **kwds):
        """Initialize the Space model.

        """
# THESE ARE THE OLD TESTS, BEFORE THE CHANGE THAT NOW HAS ELEVATION CHANGES
# ONLY APPLIED TO CORE NODES:
#        array([ 0.50005858,  0.5       ,  0.5       ,  0.5       ,  0.5       ,
#            0.5       ,  0.31524982,  0.43663631,  0.48100988,  0.5       ,
#            0.5       ,  0.43662792,  0.43661476,  0.48039544,  0.5       ,
#            0.5       ,  0.48085233,  0.48039228,  0.47769742,  0.5       ,
#            0.5       ,  0.5       ,  0.5       ,  0.5       ,  0.5       ])
#        array([ 0.02316337,  1.53606698,  2.5727653 ,  3.51126678,  4.56077707,
#            1.58157495,  0.24386828,  0.37232163,  0.42741854,  5.50969486,
#            2.54008677,  0.37232688,  0.37232168,  0.42797088,  6.52641123,
#            3.55874171,  0.42756557,  0.42797367,  0.44312812,  7.55334077,
#            4.55922478,  5.5409473 ,  6.57035008,  7.5038935 ,  8.51034357])
        #assign class variables to grid fields; create necessary fields
        self.flow_receivers = grid.at_node['flow__receiver_node']
        self.stack = grid.at_node['flow__upstream_node_order']
        self.topographic__elevation = grid.at_node['topographic__elevation']
        self.slope = grid.at_node['topographic__steepest_slope']
        self.link_to_reciever = grid.at_node['flow__link_to_receiver_node']
        self.cell_area_at_node = grid.cell_area_at_node

        try:
            self.soil__depth = grid.at_node['soil__depth']
        except KeyError:
            self.soil__depth = grid.add_zeros(
                'soil__depth', at='node', dtype=float)

        if isinstance(grid, RasterModelGrid):
            self.link_lengths = grid.length_of_d8
        else:
            self.link_lengths = grid.length_of_link

        try:
            self.bedrock__elevation = grid.at_node['bedrock__elevation']
        except KeyError:
            self.bedrock__elevation = grid.add_zeros(
                'bedrock__elevation', at='node', dtype=float)
            self.bedrock__elevation[:] = self.topographic__elevation -\
                self.soil__depth
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
        # power "m", deposition rate, and erosion rates of alluvium and rock
        self.qs_in = np.zeros(grid.number_of_nodes)
        self.Q_to_the_m = np.zeros(grid.number_of_nodes)
        self.S_to_the_n = np.zeros(grid.number_of_nodes)
        self.depo_rate = np.zeros(grid.number_of_nodes)
        self.Es = np.zeros(grid.number_of_nodes)
        self.Er = np.zeros(grid.number_of_nodes)

        # store other parameters
        self.m_sp = float(m_sp)
        self.n_sp = float(n_sp)
        self.F_f = float(F_f)
        self.phi = float(phi)
        self.H_star = float(H_star)
        self.v_s = float(v_s)
        self.dt_min = dt_min

        #K's and critical values can be floats, grid fields, or arrays
        if isinstance(K_sed, string_types):
            self.K_sed = self._grid.at_node[K_sed]
        elif isinstance(K_sed, (float, int)):  # a float
            self.K_sed = float(K_sed)
        elif len(K_sed) == self.grid.number_of_nodes:
            self.K_sed = np.array(K_sed)
        else:
            raise TypeError('Supplied type of K_sed ' +
                            'was not recognised, or array was ' +
                            'not nnodes long!')

        if isinstance(K_br, string_types):
            self.K_br = self._grid.at_node[K_br]
        elif isinstance(K_br, (float, int)):  # a float
            self.K_br = float(K_br)
        elif len(K_br) == self.grid.number_of_nodes:
            self.K_br = np.array(K_br)
        else:
            raise TypeError('Supplied type of K_br ' +
                            'was not recognised, or array was ' +
                            'not nnodes long!')

        if sp_crit_sed is not None:
            if isinstance(sp_crit_sed, string_types):
                self.sp_crit_sed = self._grid.at_node[sp_crit_sed]
            elif isinstance(sp_crit_sed, (float, int)):  # a float
                self.sp_crit_sed = float(sp_crit_sed)
            elif len(sp_crit_sed) == self.grid.number_of_nodes:
                self.sp_crit_sed = np.array(sp_crit_sed)
            else:
                raise TypeError('Supplied type of sp_crit_sed ' +
                                'was not recognised, or array was ' +
                                'not nnodes long!')

        if sp_crit_br is not None:
            if isinstance(sp_crit_br, string_types):
                self.sp_crit_br = self._grid.at_node[sp_crit_br]
            elif isinstance(sp_crit_br, (float, int)):  # a float
                self.sp_crit_br = float(sp_crit_br)
                self.sp_crit_br = np.array(sp_crit_br)
            else:
                raise TypeError('Supplied type of sp_crit_br ' +
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
            self.calc_ero_rates = self.simple_stream_power
        elif self.method == 'threshold_stream_power':
            self.calc_ero_rates = self.threshold_stream_power
        elif self.method == 'stochastic_hydrology':
            self.calc_ero_rates = self.stochastic_hydrology
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
            if self.phi >= 1.0:
                raise ValueError('Porosity phi must be < 1.0')
            self.porosity_factor = 1.0 / (1.0 - self.phi)
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
        if self.discharge_method == None:
            self.Q_to_the_m[:] = np.power(self.grid.at_node['drainage_area'], self.m_sp)
        else:
            if self.discharge_method == 'area_field':
                if self.area_field is not None:
                    if isinstance(self.area_field, string_types):
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
                    if isinstance(self.discharge_field, string_types):
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
                raise TypeError('Unrecognized discharge_method')

        #TODO: FIGURE OUT WHY WE BOTHER TO CALC AND STORE BOTH ES/ER AND
        #**_EROSION_TERM?
        self.Es[:] = self.K_sed * self.Q_to_the_m * np.power(self.slope, self.n_sp) * \
            (1.0 - np.exp(-self.soil__depth / self.H_star))
        self.Er[:] = self.K_br * self.Q_to_the_m * np.power(self.slope, self.n_sp) * \
            np.exp(-self.soil__depth / self.H_star)
        self.sed_erosion_term = self.K_sed * self.Q_to_the_m * \
            np.power(self.slope, self.n_sp)
        self.br_erosion_term = self.K_br * self.Q_to_the_m * \
            np.power(self.slope, self.n_sp)

    def threshold_stream_power(self):
        """Use stream power with entrainment/erosion thresholds.

        threshold_stream_power works the same way as simple SP but includes
        user-defined thresholds for sediment entrainment and bedrock erosion.
        """
        self.Q_to_the_m = np.zeros(len(self.grid.at_node['drainage_area']))
        if self.method == 'threshold_stream_power' and self.discharge_method == None:
            self.Q_to_the_m[:] = np.power(self.grid.at_node['drainage_area'], self.m_sp)
        elif self.method == 'threshold_stream_power' and self.discharge_method is not None:
            if self.discharge_method == 'drainage_area':
                if self.area_field is not None:
                    if isinstance(self.area_field, string_types):
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
                    if isinstance(self.discharge_field, string_types):
                        self.q[:] = self._grid.at_node[self.discharge_field]
                        self.Q_to_the_m[:] = np.power(self.q, self.m_sp)
                    elif len(self.discharge_field) == self.grid.number_of_nodes:
                        self.q[:] = np.array(self.discharge_field)
                        self.Q_to_the_m[:] = np.power(self.q, self.m_sp)
                    else:
                        raise TypeError('Supplied type of discharge_field ' +
                                'was not recognised, or array was ' +
                                'not nnodes long!')
        omega_sed = self.K_sed * self.Q_to_the_m * \
            np.power(self.slope, self.n_sp)
        omega_br = self.K_br * self.Q_to_the_m * \
            np.power(self.slope, self.n_sp)
        self.Es = (omega_sed - self.sp_crit_sed * (1 - np.exp(-omega_sed /\
            self.sp_crit_sed))) * \
            (1.0 - np.exp(-self.soil__depth / self.H_star))
        self.Er = (omega_br - self.sp_crit_br * (1 - np.exp(-omega_br /\
            self.sp_crit_br))) * \
            np.exp(-self.soil__depth / self.H_star)
        self.sed_erosion_term = omega_sed - self.sp_crit_sed * \
            (1 - np.exp(-omega_sed / self.sp_crit_sed))
        self.br_erosion_term = omega_br - self.sp_crit_br * \
            (1 - np.exp(-omega_br / self.sp_crit_br))

    def stochastic_hydrology(self):
        """Allows custom area and discharge fields, no default behavior.

        stochastic_hydrology forces the user to supply either an array or
        field name for either drainage area or discharge, and will not
        default to q=A^m.
        """
        self.Q_to_the_m = np.zeros(len(self.grid.at_node['drainage_area']))
        if self.method == 'stochastic_hydrology' and self.discharge_method == None:
            raise TypeError('Supply a discharge method to use stoc. hydro!')
        elif self.discharge_method is not None:
            if self.discharge_method == 'drainage_area':
                if self.area_field is not None:
                    if isinstance(self.area_field, string_types):
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
                    if isinstance(self.discharge_field, string_types):
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
        self.Es = self.K_sed * self.Q_to_the_m * np.power(self.slope, self.n_sp) * \
            (1.0 - np.exp(-self.soil__depth / self.H_star))
        self.Er = self.K_br * self.Q_to_the_m * np.power(self.slope, self.n_sp) * \
            np.exp(-self.soil__depth / self.H_star)
        self.sed_erosion_term = self.K_sed * self.Q_to_the_m * \
            np.power(self.slope, self.n_sp)
        self.br_erosion_term = self.K_br * self.Q_to_the_m * \
            np.power(self.slope, self.n_sp)

    def run_one_step_basic(self, dt=1.0, flooded_nodes=None, **kwds):
        """Calculate change in rock and alluvium thickness for
        a time period 'dt'.

        Parameters
        ----------
        dt : float
            Model timestep [T]
        flooded_nodes : array
            Indices of flooded nodes, passed from flow router
        """
        #Choose a method for calculating erosion:
        if self.method == 'stochastic_hydrology':
            self.stochastic_hydrology()
        elif self.method == 'simple_stream_power':
            self.simple_stream_power()
        elif self.method == 'threshold_stream_power':
            self.threshold_stream_power()
        else:
            raise ValueError('Specify an erosion method!')

        self.qs_in[:] = 0

        #iterate top to bottom through the stack, calculate qs
        # cythonized version of calculating qs_in
        calculate_qs_in(np.flipud(self.stack),
                        self.flow_receivers,
                        self.cell_area_at_node,
                        self.q,
                        self.qs,
                        self.qs_in,
                        self.Es,
                        self.Er,
                        self.v_s,
                        self.F_f,
                        self.phi)

        self.depo_rate[self.q > 0] = (self.qs[self.q > 0]
                                      * (self.v_s / self.q[self.q > 0]))

        #now, the analytical solution to soil thickness in time:
        #need to distinguish D=kqS from all other cases to save from blowup!

        flooded = np.full(self._grid.number_of_nodes, False, dtype=bool)
        flooded[flooded_nodes] = True

        #distinguish cases:
        blowup = self.depo_rate == self.K_sed * self.Q_to_the_m * self.slope

        ##first, potential blowup case:
        #positive slopes, not flooded
        self.soil__depth[(self.q > 0) & (blowup==True) & (self.slope > 0) & \
            (flooded==False)] = self.H_star * \
            np.log((self.sed_erosion_term[(self.q > 0) & (blowup==True) & \
            (self.slope > 0) & (flooded==False)] / self.H_star) * dt + \
            np.exp(self.soil__depth[(self.q > 0) & (blowup==True) & \
            (self.slope > 0) & (flooded==False)] / self.H_star))
        #positive slopes, flooded
        self.soil__depth[(self.q > 0) & (blowup==True) & (self.slope > 0) & \
            (flooded==True)] = (self.depo_rate[(self.q > 0) & \
            (blowup==True) & (self.slope > 0) & (flooded==True)] / (1 - self.phi)) * dt
        #non-positive slopes, not flooded
        self.soil__depth[(self.q > 0) & (blowup==True) & (self.slope <= 0) & \
            (flooded==False)] += (self.depo_rate[(self.q > 0) & \
            (blowup==True) & (self.slope <= 0) & (flooded==False)] / \
            (1 - self.phi)) * dt

        ##more general case:
        #positive slopes, not flooded
        self.soil__depth[(self.q > 0) & (blowup==False) & (self.slope > 0) & \
            (flooded==False)] = self.H_star * \
            np.log((1 / ((self.depo_rate[(self.q > 0) & (blowup==False) & \
            (self.slope > 0) & (flooded==False)] / (1 - self.phi)) / \
            (self.sed_erosion_term[(self.q > 0) & (blowup==False) & \
            (self.slope > 0) & (flooded==False)]) - 1)) * \
            (np.exp((self.depo_rate[(self.q > 0) & (blowup==False) & \
            (self.slope > 0) & (flooded==False)] / (1 - self.phi) - \
            (self.sed_erosion_term[(self.q > 0) & (blowup==False) & \
            (self.slope > 0) & (flooded==False)]))*(dt / self.H_star)) * \
            (((self.depo_rate[(self.q > 0) & (blowup==False) & \
            (self.slope > 0) & (flooded==False)] / (1 - self.phi) / \
            (self.sed_erosion_term[(self.q > 0) & (blowup==False) & \
            (self.slope > 0) & (flooded==False)])) - 1) * \
            np.exp(self.soil__depth[(self.q > 0) & (blowup==False) & \
            (self.slope > 0) & (flooded==False)] / self.H_star)  + 1) - 1))
        #places where slope <= 0 but not flooded:
        self.soil__depth[(self.q > 0) & (blowup==False) & (self.slope <= 0) & \
            (flooded==False)] += (self.depo_rate[(self.q > 0) & \
            (blowup==False) & (self.slope <= 0) & (flooded==False)] / \
            (1 - self.phi)) * dt
        #flooded nodes:
        self.soil__depth[(self.q > 0) & (blowup==False) & (flooded==True)] += \
            (self.depo_rate[(self.q > 0) & (blowup==False) & \
            (flooded==True)] / (1 - self.phi)) * dt

        self.bedrock__elevation[self.q > 0] += dt * \
            (-self.br_erosion_term[self.q > 0] * \
            (np.exp(-self.soil__depth[self.q > 0] / self.H_star)))

        #finally, determine topography by summing bedrock and soil
#        self.topographic__elevation[:] = self.bedrock__elevation + \
#            self.soil__depth
        cores = self._grid.core_nodes
        self.topographic__elevation[cores] = self.bedrock__elevation[cores] + \
            self.soil__depth[cores]

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
        >>> sp = Space(rg, K_sed=0.00001, K_br=0.00000000001,\
                                F_f=0.5, phi=0.1, H_star=1., v_s=0.001,\
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

    def run_with_adaptive_time_step_solver(self, dt=1.0, flooded_nodes=[],
                                           **kwds):
        """Run step with CHILD-like solver that adjusts time steps to prevent
        slope flattening.

        Examples
        --------
        >>> from landlab import RasterModelGrid
        >>> from landlab.components import FlowAccumulator
        >>> import numpy as np

        >>> rg = RasterModelGrid((3, 4))
        >>> z = rg.add_zeros('topographic__elevation', at='node')
        >>> z[:] = 0.1 * rg.x_of_node
        >>> H = rg.add_zeros('soil__depth', at='node')
        >>> H += 0.1
        >>> br = rg.add_zeros('bedrock__elevation', at='node')
        >>> br[:] = z - H

        >>> fa = FlowAccumulator(rg, flow_director='FlowDirectorSteepest')
        >>> fa.run_one_step()
        >>> sp = Space(rg, K_sed=1.0, K_br=0.1,
        ...            F_f=0.5, phi=0.0, H_star=1., v_s=1.0,
        ...            m_sp=0.5, n_sp = 1.0, sp_crit_sed=0,
        ...            sp_crit_br=0, method='simple_stream_power',
        ...            discharge_method='area_field',
        ...            area_field='drainage_area',
        ...            discharge_field=None, solver='adaptive')
        >>> sp.run_one_step(dt=10.0)

        >>> np.round(sp.Es[5:7], 4)
        array([ 0.0029,  0.0074])
        >>> np.round(sp.Er[5:7], 4)
        array([ 0.0032,  0.0085])
        >>> np.round(H[5:7], 3)
        array([ 0.088,  0.078])
        """

        # Initialize remaining_time, which records how much of the global time
        # step we have yet to use up.
        remaining_time = dt

        z = self._grid.at_node['topographic__elevation']
        br = self._grid.at_node['bedrock__elevation']
        H = self._grid.at_node['soil__depth']
        r = self.flow_receivers
        time_to_flat = np.zeros(len(z))
        time_to_zero_alluv = np.zeros(len(z))
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
                # of the dynamic timestepper, but just in case, we update here.
                new_flooded_nodes = np.where(self.slope<0)[0]
                flooded_nodes = np.asarray(np.unique(np.concatenate((
                        flooded_nodes, new_flooded_nodes))), dtype=np.int64)
            else:
                first_iteration = False

            # Calculate rates of entrainment
            self.calc_ero_rates()
            self.Es[flooded_nodes] = 0.0
            self.Er[flooded_nodes] = 0.0

            # Zero out sediment influx for new iteration
            self.qs_in[:] = 0.0

            calculate_qs_in(np.flipud(self.stack),
                            self.flow_receivers,
                            self.cell_area_at_node,
                            self.q,
                            self.qs,
                            self.qs_in,
                            self.Es,
                            self.Er,
                            self.v_s,
                            self.F_f,
                            self.phi)

            self.depo_rate[self.q > 0] = (self.qs[self.q > 0]
                                          * (self.v_s / self.q[self.q > 0]))
            # TODO handle flooded nodes in the above fn

            # Now look at upstream-downstream node pairs, and recording the
            # time it would take for each pair to flatten. Take the minimum.
            dzdt[cores] = self.depo_rate[cores] - (self.Es[cores] + self.Er[cores])
            rocdif = dzdt - dzdt[r]
            zdif = z - z[r]
            time_to_flat[:] = remaining_time

            converging = np.where(rocdif < 0.0)[0]
            time_to_flat[converging] = -(TIME_STEP_FACTOR * zdif[converging]
                                         / rocdif[converging])
            time_to_flat[np.where(zdif <= 0.0)[0]] = remaining_time

            # From this, find the maximum stable time step with regard to slope
            # evolution.
            dt_max1 = np.amin(time_to_flat)

            # Next we consider time to exhaust regolith
            time_to_zero_alluv[:] = remaining_time
            dHdt = self.porosity_factor * (self.depo_rate - self.Es)
            decreasing_H = np.where(dHdt < 0.0)[0]
            time_to_zero_alluv[decreasing_H] = - (TIME_STEP_FACTOR
                                                  * H[decreasing_H]
                                                  / dHdt[decreasing_H])

            # Now find the smallest time that would lead to near-empty alluv
            dt_max2 = np.amin(time_to_zero_alluv)

            # Take the smaller of the limits
            dt_max = min(dt_max1, dt_max2)
            if dt_max < self.dt_min:
                dt_max = self.dt_min

            # Now a vector operation: apply dzdt and dhdt to all nodes
            br[cores] -= self.Er[cores] * dt_max
            H[cores] += dHdt[cores] * dt_max
            z[cores] = br[cores] + H[cores]

            # Update remaining time and continue
            remaining_time -= dt_max

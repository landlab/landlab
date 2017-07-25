import numpy as np
from landlab import Component
from .cfuncs import calculate_qs_in

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
                 m_sp=None, n_sp=None, sp_crit=None, 
                 method=None, discharge_method=None, 
                 area_field=None, discharge_field=None, **kwds):
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
        
        >>> mg.at_node['topographic__elevation'] # doctest: +NORMALIZE_WHITESPACE
        array([-0.47709402,  1.03606698,  2.0727653 ,  3.01126678,  4.06077707,
            1.08157495, -0.0799798 , -0.06459322, -0.05380581,  5.00969486,
            2.04008677, -0.06457996, -0.06457219, -0.05266169,  6.02641123,
            3.05874171, -0.05350698, -0.05265586, -0.03498794,  7.05334077,
            4.05922478,  5.0409473 ,  6.07035008,  7.0038935 ,  8.01034357])        
        """
                #assign class variables to grid fields; create necessary fields
        self.flow_receivers = grid.at_node['flow__receiver_node']
        self.stack = grid.at_node['flow__upstream_node_order']
        self.elev = grid.at_node['topographic__elevation']
        self.slope = grid.at_node['topographic__steepest_slope']
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
        
        #store other constants
        self.m_sp = float(m_sp)
        self.n_sp = float(n_sp)
        self.phi = float(phi)
        self.v_s = float(v_s)
        
        #K's and critical values can be floats, grid fields, or arrays
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
            self.simple_stream_power()
        elif self.method == 'threshold_stream_power':
            self.threshold_stream_power()
        elif self.method == 'stochastic_hydrology':
            self.stochastic_hydrology()
        else:
            raise ValueError('Specify erosion method (simple stream power,\
                            threshold stream power, or stochastic hydrology)!')
    #three choices for erosion methods:
    def simple_stream_power(self):
        """Use non-threshold stream power.
        
        simple_stream_power uses no entrainment or erosion thresholds,
        and uses either q=A^m or q=Q^m depending on discharge method. If
        discharge method is None, default is q=A^m.
        """
        self.Q_to_the_m = np.zeros(len(self.grid.at_node['drainage_area']))
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
        self.erosion_term = self.K * self.Q_to_the_m * \
            np.power(self.slope, self.n_sp)
        self.qs_in = np.zeros(self.grid.number_of_nodes) 
            
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
        omega = self.K * self.Q_to_the_m * np.power(self.slope, self.n_sp)
        self.erosion_term = omega - self.sp_crit * \
            (1 - np.exp(-omega / self.sp_crit))
        self.qs_in = np.zeros(self.grid.number_of_nodes) 
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
        self.erosion_term = self.K * self.Q_to_the_m * \
            np.power(self.slope, self.n_sp)
        self.qs_in = np.zeros(self.grid.number_of_nodes) 

    def run_one_step(self, dt=1.0, flooded_nodes=None, **kwds):
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
            
        self.qs_in[:] = 0# np.zeros(self.grid.number_of_nodes)            
        #iterate top to bottom through the stack, calculate qs
        # cythonized version of calculating qs_in
        calculate_qs_in(np.flipud(self.stack),
                        self.flow_receivers,
                        self.grid.node_spacing,
                        self.q,
                        self.qs,
                        self.qs_in,
                        self.erosion_term,
                        self.v_s)
    
        deposition_pertime = np.zeros(self.grid.number_of_nodes)
        deposition_pertime[self.q > 0] = (self.qs[self.q > 0] * \
                                         (self.v_s / self.q[self.q > 0]))

        #now, the analytical solution to soil thickness in time:
        #need to distinguish D=kqS from all other cases to save from blowup!
        
        flooded = np.full(self._grid.number_of_nodes, False, dtype=bool)
        flooded[flooded_nodes] = True        
        
        #topo elev is old elev + deposition - erosion        
        self.elev += ((deposition_pertime - self.erosion_term) * dt)        
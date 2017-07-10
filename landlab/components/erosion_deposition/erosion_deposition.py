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
        """Initialize the HybridAlluvium model.
        
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
                raise TypeError('Supplied type of sp_crit_sed ' +
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
        self.erosion_term = omega - self.sp_crit_sed * \
            (1 - np.exp(-omega / self.sp_crit_sed))
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
                        self.Es,
                        self.Er,
                        self.v_s,
                        self.F_f)
    
        deposition_pertime = np.zeros(self.grid.number_of_nodes)
        deposition_pertime[self.q > 0] = (self.qs[self.q > 0] * \
                                         (self.v_s / self.q[self.q > 0]))

        #now, the analytical solution to soil thickness in time:
        #need to distinguish D=kqS from all other cases to save from blowup!
        
        flooded = np.full(self._grid.number_of_nodes, False, dtype=bool)
        flooded[flooded_nodes] = True        
        
        #distinguish cases:
        blowup = deposition_pertime == self.K_sed * self.Q_to_the_m * self.slope

        ##first, potential blowup case:
        #positive slopes, not flooded
        self.soil__depth[(self.q > 0) & (blowup==True) & (self.slope > 0) & \
            (flooded==False)] = self.H_star * \
            np.log((self.erosion_term[(self.q > 0) & (blowup==True) & \
            (self.slope > 0) & (flooded==False)] / self.H_star) * dt + \
            np.exp(self.soil__depth[(self.q > 0) & (blowup==True) & \
            (self.slope > 0) & (flooded==False)] / self.H_star))
        #positive slopes, flooded
        self.soil__depth[(self.q > 0) & (blowup==True) & (self.slope > 0) & \
            (flooded==True)] = (deposition_pertime[(self.q > 0) & \
            (blowup==True) & (self.slope > 0) & (flooded==True)] / (1 - self.phi)) * dt   
        #non-positive slopes, not flooded
        self.soil__depth[(self.q > 0) & (blowup==True) & (self.slope <= 0) & \
            (flooded==False)] += (deposition_pertime[(self.q > 0) & \
            (blowup==True) & (self.slope <= 0) & (flooded==False)] / \
            (1 - self.phi)) * dt    
        
        ##more general case:
        #positive slopes, not flooded
        self.soil__depth[(self.q > 0) & (blowup==False) & (self.slope > 0) & \
            (flooded==False)] = self.H_star * \
            np.log((1 / ((deposition_pertime[(self.q > 0) & (blowup==False) & \
            (self.slope > 0) & (flooded==False)] / (1 - self.phi)) / \
            (self.erosion_term[(self.q > 0) & (blowup==False) & \
            (self.slope > 0) & (flooded==False)]) - 1)) * \
            (np.exp((deposition_pertime[(self.q > 0) & (blowup==False) & \
            (self.slope > 0) & (flooded==False)] / (1 - self.phi) - \
            (self.erosion_term[(self.q > 0) & (blowup==False) & \
            (self.slope > 0) & (flooded==False)]))*(dt / self.H_star)) * \
            (((deposition_pertime[(self.q > 0) & (blowup==False) & \
            (self.slope > 0) & (flooded==False)] / (1 - self.phi) / \
            (self.erosion_term[(self.q > 0) & (blowup==False) & \
            (self.slope > 0) & (flooded==False)])) - 1) * \
            np.exp(self.soil__depth[(self.q > 0) & (blowup==False) & \
            (self.slope > 0) & (flooded==False)] / self.H_star)  + 1) - 1))
        #places where slope <= 0 but not flooded:
        self.soil__depth[(self.q > 0) & (blowup==False) & (self.slope <= 0) & \
            (flooded==False)] += (deposition_pertime[(self.q > 0) & \
            (blowup==False) & (self.slope <= 0) & (flooded==False)] / \
            (1 - self.phi)) * dt     
        #flooded nodes:        
        self.soil__depth[(self.q > 0) & (blowup==False) & (flooded==True)] += \
            (deposition_pertime[(self.q > 0) & (blowup==False) & \
            (flooded==True)] / (1 - self.phi)) * dt     
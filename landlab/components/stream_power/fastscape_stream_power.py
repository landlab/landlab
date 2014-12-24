#! /usr/env/python

"""
This module attempts to "component-ify" GT's Fastscape stream power erosion.
Created DEJH, March 2014.
"""

import numpy
from landlab import ModelParameterDictionary
from landlab.core.model_parameter_dictionary import MissingKeyError, ParameterValueError
from landlab.field.scalar_data_fields import FieldError
from scipy import weave
#from scipy.weave.build_tools import CompileError
UNDEFINED_INDEX = numpy.iinfo(numpy.int32).max

class SPEroder(object):
    '''
    This class uses the Braun-Willett Fastscape approach to calculate the amount
    of erosion at each node in a grid, following a stream power framework.
    
    On initialization, it takes *grid*, a reference to a ModelGrid, and
    *input_stream*, a string giving the filename (and optionally, path) of the
    required input file.
    
    It needs to be supplied with the key variables:
    
        *K_sp*
        
        *m_sp*
    
    ...which it will draw from the supplied input file. *n_sp* has to be 1 for 
    the BW algorithm to work. If you want n!=1., try calling the explicit
    "stream_power" component.
    
    If you want to supply a spatial variation in K, set K_sp to the string
    'array', and pass a field name or array to the erode method's K_if_used
    argument.
    
    *dt*, *rainfall_intensity*, and *value_field* are optional variables.
    
    *dt* is a fixed timestep, and *rainfall_intensity* is a parameter which 
    modulates K_sp (by a product, r_i**m_sp) to reflect the direct influence of
    rainfall intensity on erosivity. *value_field* is a string giving the name
    of the field containing the elevation data in the grid. It defaults to
    'topographic_elevation' if not supplied.
    
    This module assumes you have already run 
    :func:`landlab.components.flow_routing.route_flow_dn.FlowRouter.route_flow`
    in the same timestep. It looks for 'upstream_ID_order', 
    'links_to_flow_receiver', 'drainage_area', 'flow_receiver', and
    'topographic_elevation' at the nodes in the grid. 'drainage_area' should
    be in area upstream, not volume (i.e., set runoff_rate=1.0 when calling
    FlowRouter.route_flow).
    
    If dt is not supplied, you must call gear_timestep(dt_in, rain_intensity_in)
    each iteration to set these variables on-the-fly (rainfall_intensity will be
    overridden if supplied in the input file).
    If dt is supplied but rainfall_intensity is not, the module will assume you
    mean r_i = 1.
    
    The primary method of this class is :func:`erode`.
    '''
    
    def __init__(self, grid, input_stream):
        self.grid = grid
        inputs = ModelParameterDictionary(input_stream)
        
        #User sets:
        try:
            self.K = inputs.read_float('K_sp')
        except ParameterValueError:
            self.use_K = True
        else:
            self.use_K = False
            
        self.m = inputs.read_float('m_sp')
        try:
            self.n = inputs.read_float('n_sp')
        except:
            self.n = 1.
        try:
            self.dt = inputs.read_float('dt')
        except: #if dt isn't supplied, it must be set by another module, so look in the grid
            print 'Set dynamic timestep from the grid. You must call gear_timestep() to set dt each iteration.'
        try:
            self.r_i = inputs.read_float('rainfall_intensity')
        except:
            self.r_i = 1.
        try:
            self.value_field = inputs.read_str('value_field')
        except:
            self.value_field = 'topographic_elevation'
            
        #make storage variables
        self.A_to_the_m = grid.create_node_array_zeros()
        self.alpha = grid.empty(centering='node')
        
        self.grid.node_diagonal_links() #calculates the number of diagonal links
        
        if self.n != 1.:
            raise ValueError('The Braun Willett stream power algorithm requires n==1. at the moment, sorry...')
        
        self.weave_flag = grid.weave_flag

    def gear_timestep(self, dt_in, rainfall_intensity_in=None):
        self.dt = dt_in
        if rainfall_intensity_in is not None:
            self.r_i = rainfall_intensity_in
        return self.dt, self.r_i

    def erode(self, grid_in, K_if_used=None):
        """
        This method implements the stream power erosion, following the Braun-
        Willett (2013) implicit Fastscape algorithm. This should allow it to
        be stable against larger timesteps than an explicit stream power scheme.
        
        The method takes *grid*, a reference to the model grid.
        Set 'K_if_used' as a field name or nnodes-long array if you set K_sp as
        'array' during initialization.
        
        It returns the grid, in which it will have modified the value of 
        *value_field*, as specified in component initialization.
        """
        
        #self.grid = grid_in #the variables must be stored internally to the grid, in fields
        upstream_order_IDs = self.grid['node']['upstream_ID_order']
        #ordered_receivers = self.grid['node']['flow_receiver'][upstream_order_IDs]  #"j" in GT's sketch
        #nonboundaries = numpy.not_equal(upstream_order_IDs, ordered_receivers)
        z = self.grid['node'][self.value_field]
        #interior_nodes = numpy.greater_equal(self.grid['node']['links_to_flow_receiver'], -1)
        #interior_nodes = (self.grid['node']['links_to_flow_receiver'][upstream_order_IDs])[nonboundaries]
        #flow_link_lengths = self.grid.link_length[interior_nodes]
        ##defined_flow_receivers = numpy.greater_equal(self.grid['node']['links_to_flow_receiver'],-1)
        defined_flow_receivers = numpy.not_equal(self.grid['node']['links_to_flow_receiver'],UNDEFINED_INDEX)
        #flow_link_lengths = numpy.zeros_like(self.alpha)
        flow_link_lengths = self.grid.link_length[self.grid['node']['links_to_flow_receiver'][defined_flow_receivers]]
        
        if K_if_used!=None:
            assert self.use_K, "An array of erodabilities was provided, but you didn't set K_sp to 'array' in your input file! Aborting..."
            try:
                self.K = self.grid.at_node[K_if_used][defined_flow_receivers]
            except TypeError:
                self.K = K_if_used[defined_flow_receivers]
        
        #regular_links = numpy.less(self.grid['node']['links_to_flow_receiver'][defined_flow_receivers],self.grid.number_of_links)
        #flow_link_lengths[defined_flow_receivers][regular_links] = self.grid.link_length[(self.grid['node']['links_to_flow_receiver'])[defined_flow_receivers][regular_links]]
        #diagonal_links = numpy.logical_not(regular_links)
        #flow_link_lengths[defined_flow_receivers][diagonal_links] = numpy.sqrt(self.grid.node_spacing*self.grid.node_spacing)
        numpy.power(self.grid['node']['drainage_area'], self.m, out=self.A_to_the_m)
        #self.alpha[nonboundaries] = self.K * self.dt * self.A_to_the_m[nonboundaries] / flow_link_lengths
        self.alpha[defined_flow_receivers] = self.r_i**self.m * self.K * self.dt * self.A_to_the_m[defined_flow_receivers] / flow_link_lengths

        flow_receivers = self.grid['node']['flow_receiver']
        n_nodes = upstream_order_IDs.size
        alpha = self.alpha
        if self.weave_flag:
            code = """
                int current_node;
                int j;
                for (int i = 0; i < n_nodes; i++) {
                    current_node = upstream_order_IDs[i];
                    j = flow_receivers[current_node];
                    if (current_node != j) {
                        z[current_node] = (z[current_node] + alpha[current_node]*z[j])/(1.0+alpha[current_node]);
                    }
                }
            """
            weave.inline(code, ['n_nodes', 'upstream_order_IDs', 'flow_receivers', 'z', 'alpha'])
        else:
            for i in upstream_order_IDs:
                j = flow_receivers[i]
                if i != j:
                    z[i] = (z[i] + alpha[i]*z[j])/(1.0+alpha[i])
        
        #self.grid['node'][self.value_field] = z
        
        return self.grid


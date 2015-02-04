#! /usr/env/python

"""
This module attempts to "component-ify" GT's Fastscape stream power erosion.
Created DEJH, March 2014.
"""

import numpy
from landlab import ModelParameterDictionary
from landlab import Component
from landlab.core.model_parameter_dictionary import MissingKeyError, ParameterValueError
from landlab.field.scalar_data_fields import FieldError
from scipy import weave
from scipy.optimize import newton, fsolve
#from scipy.weave.build_tools import CompileError
UNDEFINED_INDEX = numpy.iinfo(numpy.int32).max

class SPEroder(Component):
    '''
    This class uses the Braun-Willett Fastscape approach to calculate the amount
    of erosion at each node in a grid, following a stream power framework.
    
    On initialization, it takes *grid*, a reference to a ModelGrid, and
    *input_stream*, a string giving the filename (and optionally, path) of the
    required input file.
    
    It needs to be supplied with the key variables:
    
        *K_sp*
        *m_sp*
        *n_sp*
    
    ...which it will draw from the supplied input file. *n_sp*  can be any 
    value ~ 0.5<n_sp<4., but note that performance will be EXTREMELY degraded
    if n<1.
    
    If you want to supply a spatial variation in K, set K_sp to the string
    'array', and pass a field name or array to the erode method's K_if_used
    argument.
    
    *rainfall_intensity* is an optional variable; a parameter which 
    modulates K_sp (by a product, r_i**m_sp) to reflect the direct influence of
    rainfall intensity on erosivity.
    
    This module assumes you have already run 
    :func:`landlab.components.flow_routing.route_flow_dn.FlowRouter.route_flow`
    in the same timestep. It looks for 'upstream_ID_order', 
    'links_to_flow_receiver', 'drainage_area', 'flow_receiver', and
    'topographic_elevation' at the nodes in the grid, or needs to be supplied
    with them directly when running. 'drainage_area' should
    be in area upstream, not volume (i.e., set runoff_rate=1.0 when calling
    FlowRouter.route_flow).
    
    The primary method of this class is :func:`erode`.
    '''

    _name = 'fastscape_stream_power'
    
    _input_var_names = set([
        'topographic_elevation',
        'drainage_area', 
        'upstream_ID_order', 
        'links_to_flow_receiver', 
        'flow_receiver'
    ])
    
    _output_var_names = set([
        'topographic_elevation'
    ])
    
    _var_units = {
        'topographic_elevation' : 'm',
        'drainage_area' : 'm**2',
        'upstream_ID_order' : 'None',
        'links_to_flow_receiver' : 'None',
        'flow_receiver' : 'None'
    }
    
    _var_mapping = {
        'topographic_elevation' : 'Node',
        'drainage_area' : 'Node',
        'upstream_ID_order' : 'Node',
        'links_to_flow_receiver' : 'Node',
        'flow_receiver' : 'Node'
    }
    
    _var_defs = {
        'topographic_elevation' : 'The land surface elevation',
        'drainage_area' : 'The upstream accumulated drainage area at each node',
        'upstream_ID_order' : 'The node ID of each node in the grid, arranged in order from outlet to drainage network tips by distance from outlet',
        'links_to_flow_receiver' : 'The link ID connecting each node to its downstream neighbor',
        'flow_receiver' : 'The node ID of the node downstream for each node'
    }
    
    def __init__(self, grid, input_stream):
        
        self._grid = grid
        inputs = ModelParameterDictionary(input_stream)
        
        #perform a check that none of the standard names have been overwritten in the params file:
        for current_var_names in (self._input_var_names, self._output_var_names):
            for io_name in current_var_names:
                try:
                    updated_name = inputs.read_string(io_name)
                except:
                    pass
                else:
                    del current_var_names[io_name]
                    current_var_names.add(updated_name)
                    for this_dict in (self._var_units, self._var_mapping, self._var_defs):
                        this_value = this_dict[io_name]
                        del this_dict[io_name]
                        this_dict[updated_name] = this_value
                    
        
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
        self.alpha_by_flow_link_lengthtothenless1 = numpy.empty_like(self.alpha)
        
        self._grid.node_diagonal_links() #calculates the number of diagonal links
        
        if self.n != 1.:
            #raise ValueError('The Braun Willett stream power algorithm requires n==1. at the moment, sorry...')
            self.nonlinear_flag = True
            if self.n<1.:
                print "***WARNING: With n<1 performance of the Fastscape algorithm is slow!***"
        else:
            self.nonlinear_flag = False
        
        self.weave_flag = grid.weave_flag
        
        def func_for_newton(x, last_step_elev, receiver_elev, alpha_by_flow_link_lengthtothenless1, n):
            y = x - last_step_elev + alpha_by_flow_link_lengthtothenless1*(x-receiver_elev)**n
            return y
        
        def func_for_newton_diff(x, last_step_elev, receiver_elev, alpha_by_flow_link_lengthtothenless1, n):
            y = 1. + n*alpha_by_flow_link_lengthtothenless1*(x-receiver_elev)**(n-1.)
            return y
        
        self.func_for_newton = func_for_newton
        self.func_for_newton_diff = func_for_newton_diff
        
        
    def erode(self, dt, 
              topographic_elevation='topographic_elevation',
              drainage_area='drainage_area',
              upstream_ID_order='upstream_ID_order',
              links_to_flow_receiver='links_to_flow_receiver',
              flow_receiver='flow_receiver',
              K_if_used=None, rainfall_intensity_if_used=None):
        """
        This method implements the stream power erosion, following the Braun-
        Willett (2013) implicit Fastscape algorithm. This should allow it to
        be stable against larger timesteps than an explicit stream power scheme.
        
        The method takes *dt*, the timestep.
        
        It needs access to *topographic_elevation*, *drainage_area*,
        *upstream_ID_order*, *links_to_flow_receiver*, & *flow_receiver*. These
        can be supplied either through the grid fields as strings 
        (recommended), or can be supplied direct as numpy arrays. The default
        strings reflect Landlab standard naming conventions, and the latter
        three will be set to these names by route_flow_dn.FlowRouter by
        default.
        
        Set 'K_if_used' as a field name or nnodes-long array if you set K_sp as
        'array' during initialization.
        Set 'rainfall_intensity_if_used' as a field name or nnodes-long array
        if you want spatially variable rainfall (leave as None if you specified
        a float in the input file or just want stream power determined by A).
        
        This method updates the topographic elevation in the grid field (if 
        supplied), and also returns a length-1 tuple containing the array of 
        topographic elevations.
        """
        
        self.dt=dt
        try:
            upstream_order_IDs = self._grid['node'][upstream_ID_order]
        except TypeError:
            upstream_order_IDs = upstream_ID_order
        try:
            z = self._grid['node'][topographic_elevation]
        except TypeError:
            z = topographic_elevation
        try:
            A = self._grid.at_node[drainage_area]
        except TypeError:
            A = drainage_area
        try:
            links_to_flow_receiver_in = self._grid['node'][links_to_flow_receiver]
        except:
            links_to_flow_receiver_in = links_to_flow_receiver
        defined_flow_receivers = numpy.not_equal(links_to_flow_receiver_in,UNDEFINED_INDEX)
        #flow_link_lengths = numpy.zeros_like(self.alpha)
        flow_link_lengths = self._grid.link_length[links_to_flow_receiver_in[defined_flow_receivers]]
        
        if K_if_used!=None:
            assert self.use_K, "An array of erodabilities was provided, but you didn't set K_sp to 'array' in your input file! Aborting..."
            try:
                self.K = self._grid.at_node[K_if_used][defined_flow_receivers]
            except TypeError:
                self.K = K_if_used[defined_flow_receivers]
        
        if rainfall_intensity_if_used!=None:
            try:
                self.r_i = self._grid.at_node[rainfall_intensity_if_used][defined_flow_receivers]
            except TypeError:
                self.r_i = rainfall_intensity_if_used[defined_flow_receivers]
            
        
        #regular_links = numpy.less(self._grid['node']['links_to_flow_receiver'][defined_flow_receivers],self._grid.number_of_links)
        #flow_link_lengths[defined_flow_receivers][regular_links] = self._grid.link_length[(self._grid['node']['links_to_flow_receiver'])[defined_flow_receivers][regular_links]]
        #diagonal_links = numpy.logical_not(regular_links)
        #flow_link_lengths[defined_flow_receivers][diagonal_links] = numpy.sqrt(self._grid.node_spacing*self._grid.node_spacing)
        numpy.power(A, self.m, out=self.A_to_the_m)
        #self.alpha[nonboundaries] = self.K * self.dt * self.A_to_the_m[nonboundaries] / flow_link_lengths
        self.alpha[defined_flow_receivers] = self.r_i**self.m * self.K * self.dt * self.A_to_the_m[defined_flow_receivers] / flow_link_lengths

        try:        
            flow_receivers = self._grid['node'][flow_receiver]
        except TypeError:
            flow_receivers = flow_receiver
        n_nodes = upstream_order_IDs.size
        alpha = self.alpha
        if self.nonlinear_flag==False: #n==1
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
        else: #general, nonlinear n case
            self.alpha_by_flow_link_lengthtothenless1[defined_flow_receivers] = alpha[defined_flow_receivers]/flow_link_lengths**(self.n-1.)
            alpha_by_flow_link_lengthtothenless1 = self.alpha_by_flow_link_lengthtothenless1
            n = float(self.n)
            if self.weave_flag:
                if n<1.:
                    #this is SLOOOOOOOOOOOW...
                    for i in upstream_order_IDs:
                        j = flow_receivers[i]
                        func_for_newton = self.func_for_newton
                        func_for_newton_diff = self.func_for_newton_diff
                        if i != j:
                            z[i] = fsolve(func_for_newton, z[i], args=(z[i], z[j], alpha_by_flow_link_lengthtothenless1[i], n))
                else:
                    code = """
                        int current_node;
                        int j;
                        double current_z;
                        double previous_z;
                        double elev_diff;
                        double elev_diff_tothenless1;
                        for (int i = 0; i < n_nodes; i++) {
                            current_node = upstream_order_IDs[i];
                            j = flow_receivers[current_node];
                            previous_z = z[current_node];
                            if (current_node != j) {
                                while (1) {
                                    elev_diff = previous_z-z[j];
                                    elev_diff_tothenless1 = pow(elev_diff, n-1.); //this isn't defined if in some iterations the elev_diff goes -ve
                                    current_z = previous_z - (previous_z - z[current_node] + alpha_by_flow_link_lengthtothenless1[current_node]*elev_diff_tothenless1*elev_diff)/(1.+n*alpha_by_flow_link_lengthtothenless1[current_node]*elev_diff_tothenless1);
                                    if (abs((current_z - previous_z)/current_z) < 1.48e-08) break;
                                    previous_z = current_z;
                                }
                                z[current_node] = current_z;
                            }
                        }
                    """
                    weave.inline(code, ['n_nodes', 'upstream_order_IDs', 'flow_receivers', 'z', 'alpha_by_flow_link_lengthtothenless1', 'n'], headers=["<math.h>"])
            else:
                for i in upstream_order_IDs:
                    j = flow_receivers[i]
                    func_for_newton = self.func_for_newton
                    func_for_newton_diff = self.func_for_newton_diff
                    if i != j:
                        if n>=1.:
                            z[i] = newton(func_for_newton, z[i], fprime=func_for_newton_diff, args=(z[i], z[j], alpha_by_flow_link_lengthtothenless1[i], n), maxiter=10)
                        else:
                            z[i] = fsolve(func_for_newton, z[i], args=(z[i], z[j], alpha_by_flow_link_lengthtothenless1[i], n))
        #self._grid['node'][self.value_field] = z
        
        return (z,)


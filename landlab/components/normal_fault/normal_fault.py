#!/usr/bin/env python
"""Rock uplift along a normal fault.

Landlab component that implements rock uplift by a normal fault. This component
does not make any attempt to advect topography laterally.

 """

import numpy as np
from landlab import Component
from landlab.utils.decorators import use_field_name_or_array


TWO_PI = 2.0*np.pi

@use_field_name_or_array('node')
def _return_surface(grid, surface):
    """
    Private function to return the surface modfy with the normal fault.

    This function exists to take advantange of the 'use_field_name_or_array
    decorator which permits providing the surface as a field name or array.
    """
    return surface

class NormalFault(Component):
    """NormalFault implements relative rock motion due to a normal fault.
    
    The fault can have an arbitrary trace given by two points (x1, y1) and 
    (x2, y2) in the `fault_trace_dict` input parameter. 
    
    This NormalFault component permits a fault activity onset and end time such
    that the uplift-time pattern implemented is a piecewise function with up to
    three segments. 
    
    """



    _name = 'NormalFault'

    #_cite_as = """ """

    _input_var_names = (
        'topographic__elevation',
    )

    _output_var_names = (
        'topographic__elevation',
    )

    _var_units = {
        'topographic__elevation': 'm',
    }

    _var_mapping = {
        'topographic__elevation': 'node',
    }

    _var_doc = {
        'topographic__elevation': 'elevation of the ground surface'
    }

    def __init__(self, grid, params):
        """
        Instantiation of a NormalFault. 
        
        Parameters
        --------  
        grid : ModelGrid
        faulted_surface : str or ndarray of shape `(n_nodes, )`
            Surface that 
            Default value is `topographic__elevation`
        uplift_start_time : float
            Text
        uplift_end_time : float
        pre_onset_throw_rate : float
        active_throw_rate : float
        post_stabilization_throw_rate : float
        fault_dip_angle : float
            Dip angle of the fault in 
            
        fault_trace_dict : dictionary
            = {x1:
                            y1:
                            x2:
                            y2:}
        include_boundaries : boolean
            Default value is False
        
        
        
         Examples
         --------
        
         Create a grid on which we will run the normal fault.
        
         >>> from landlab import RasterModelGrid
         >>> from landlab.components.flexure import Flexure
         >>> grid = RasterModelGrid((5, 4), spacing=(1.e4, 1.e4))
        
         Check the fields that are used as input to the flexure component.
        
         >>> Flexure.input_var_names # doctest: +NORMALIZE_WHITESPACE
         ('lithosphere__overlying_pressure_increment',)
        
         Check the units for the fields.
        
         >>> Flexure.var_units('lithosphere__overlying_pressure_increment')
         'Pa'
        
         If you are not sure about one of the input or output variables, you can
         get help for specific variables.
        
         >>> Flexure.var_help('lithosphere__overlying_pressure_increment')
         name: lithosphere__overlying_pressure_increment
         description:
           Applied pressure to the lithosphere over a time step
         units: Pa
         at: node
         intent: in
        
         >>> flex = Flexure(grid)
        
         In creating the component, a field (initialized with zeros) was added to the
         grid. Reset the interior nodes for the loading.
        
         >>> dh = grid.at_node['lithosphere__overlying_pressure_increment']
         >>> dh = dh.reshape(grid.shape)
         >>> dh[1:-1, 1:-1] = flex.gamma_mantle
        
         >>> flex.update()
        
         >>> flex.output_var_names
         ('lithosphere_surface__elevation_increment',)
         >>> flex.grid.at_node['lithosphere_surface__elevation_increment']
         ...     # doctest: +NORMALIZE_WHITESPACE
         array([ 0., 0., 0., 0.,
                 0., 1., 1., 0.,
                 0., 1., 1., 0.,
                 0., 1., 1., 0.,
                 0., 0., 0., 0.])
        """
        super(NormalFault, self).__init__(grid)

        self._grid = grid

        surface = params.get('faulted_surface', 'topographic__elevation')
        self.z = _return_surface(grid, surface)

        self.throw = params['active_throw_rate']
        self.fault_dip = params['fault_dip_angle']

        self.active_uplift_rate = self.throw * np.sin(self.fault_dip)

        self.start = params['uplift_start_time']
        self.stop = params.get('uplift_end_time', params['run_duration'])

        self.pre_onset_uplift_rate = params.get('pre_onset_throw_rate', 0.0)* np.sin(self.fault_dip)
        self.post_stabilization_uplift_rate = params.get('post_stabilization_throw_rate', 0.0) * np.sin(self.fault_dip)


        self.include_boundaries = params.get('include_boundaries', True)

        self.current_time = 0.0

        self.fault_trace_dict = params['fault_trace_dict']
        dx = self.fault_trace_dict['x2'] - self.fault_trace_dict['x1']
        dy = self.fault_trace_dict['y2'] - self.fault_trace_dict['y1']
        self.fault_azimuth = np.mod(np.arctan2(dy, dx), TWO_PI)
        self.fault_anti_azimuth = self.fault_azimuth + np.pi

        if dx == 0:
            self.dy_over_dx = 0.0
            self.fault_trace_y_intercept = 0.0
            self.fault_trace_x_intercept = self.fault_trace_dict['x2']
        else:
            self.dy_over_dx = dy/dx
            self.fault_trace_y_intercept = self.fault_trace_dict['y1'] - (self.dy_over_dx * self.fault_trace_dict['x1'])
            self.fault_trace_x_intercept = 0.0

        # select those nodes that are on the correct side of the node
        if self.include_boundaries:
            potential_nodes = np.arange(self._grid.size('node'))
        else:
            potential_nodes = self._grid.core_nodes

        dx_pn = (self._grid.x_of_node[potential_nodes] - self.fault_trace_x_intercept)
        dy_pn = (self._grid.y_of_node[potential_nodes] - self.fault_trace_y_intercept)

        potential_angles = np.mod(np.arctan2(dy_pn, dx_pn), TWO_PI)

        if self.fault_anti_azimuth <= TWO_PI:
            faulted_node_ids = potential_nodes[((potential_angles>self.fault_azimuth) &
                                                  (potential_angles <= (self.fault_anti_azimuth)))]
        else:
            faulted_node_ids = potential_nodes[((potential_angles>self.fault_azimuth) |
                                                  (potential_angles <= np.mod(self.fault_anti_azimuth, TWO_PI)))]

        # create an array to store faulted neighbors
        self.faulted_nodes = np.zeros(self._grid.size('node'), dtype=bool)
        self.faulted_nodes[faulted_node_ids] = True

    def run_one_step(self, dt):
        """Run one step method for NormalFaultHandler."""
        # If before start time, uplift and the pre-onset uplift rate
        if self.current_time < self.start:
            self.z[self.faulted_nodes] += self.pre_onset_uplift_rate * dt

        # if between start and stop, uplift at the active uplift rate
        elif self.current_time >= self.start and self.current_time < self.stop:
            self.z[self.faulted_nodes] += self.active_uplift_rate * dt

        # if after the stop time, use the post_stabilization rate.
        elif self.current_time >= self.stop:
            self.z[self.faulted_nodes] += self.post_stabilization_uplift_rate * dt

        if self.include_boundaries:
            # set faulted boundaries to average of faulted neighbors
            pass
        self.current_time += dt

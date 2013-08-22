#! /usr/env/python
"""

Component that models 2D diffusion using an explicit finite-volume method.

Created July 2013 GT
Last updated August 2013 GT

"""

from landlab import ModelParameterDictionary
from landlab import create_and_initialize_grid

_ALPHA = 0.1   # time-step stability factor

class DiffusionComponent():
    
    def __init__(self, grid=None, current_time=0.):
        
        self.grid = grid
        self.current_time = current_time
        
    def initialize(self, input_stream):
        
        # Create a ModelParameterDictionary for the inputs
        if type(input_stream)==ModelParameterDictionary:
            inputs = input_stream
        else:
            inputs = ModelParameterDictionary(input_stream)
        
        # Read input/configuration parameters
        self.kd = inputs.get('DIFMOD_KD', ptype=float)
        self.uplift_rate = inputs.get('DIFMOD_UPLIFT_RATE', 0., ptype=float)
        
        # Create grid if one doesn't already exist
        if self.grid==None:
            self.grid = create_and_initialize_grid(input_stream)
            
        # Set internal time step
        # ..todo:
        #   implement mechanism to compute time-steps dynamically if grid is 
        #   adaptive/changing
        dx = self.grid.get_minimum_active_link_length()  # smallest active link length
        self.dt = _ALPHA*dx*dx/self.kd  # CFL condition
        
        # Get a list of interior cells
        self.interior_cells = self.grid.get_active_cell_node_ids()
        
        # Create data arrays for variables that won't (?) be shared with other
        # components
        self.g = self.grid.create_active_link_dvector()  # surface gradients
        self.qs = self.grid.create_active_link_dvector()  # unit sediment flux
        self.dqds = self.grid.create_node_dvector()  # sed flux derivative
        
        
    def run_one_step(self, z, delt):
        
        # Take the smaller of delt or built-in time-step size self.dt
        dt = min(self.dt, delt)
        
        # Calculate the gradients and sediment fluxes
        self.g = self.grid.calculate_gradients_at_active_links(z)
        self.qs = -self.kd*self.g
        
        # Calculate the net deposition/erosion rate at each node
        self.dqsds = self.grid.calculate_flux_divergence_at_nodes(self.qs)
        
        # Calculate the total rate of elevation change
        dzdt = self.uplift_rate - self.dqsds
        
        # Update the elevations
        z[self.interior_cells] = z[self.interior_cells] \
                                 + dzdt[self.interior_cells] * dt
                                 
        # Update current time and return it
        self.current_time += dt
        return self.current_time
        
        
    def run_until(self, t, z):
        
        while self.current_time < t:
            remaining_time = t - self.current_time
            self.run_one_step(z, remaining_time)
            
        
    def get_time_step(self):
        """
        Returns time-step size.
        """
        return self.dt
        

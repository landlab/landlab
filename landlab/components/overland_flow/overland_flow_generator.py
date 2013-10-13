#! /usr/env/python
"""
Overland flow component.  2D numerical model of shallow-water flow over 
topography using the Bates et al. (2010) algorithm for storage-cell inundation 
modeling.

Most of code written by GT.

NG added a shear stress component.
"""

import landlab
from landlab import ModelParameterDictionary
import numpy as np
import pylab

class OverlandFlow(object):
    
    def __init__(self, input_stream, grid=None, current_time=0.):
        
        self.grid = grid
        #create and initial grid if one doesn't already exist
        if self.grid==None:
            self.grid = create_and_initialize_grid(input_stream)
        self.grid = grid
        self.current_time = current_time
        self.initialize(grid, input_stream)
        
    def initialize(self, grid, input_stream):
        
        # Create a ModelParameterDictionary for the inputs
        if type(input_stream)==ModelParameterDictionary:
            inputs = input_stream
        else:
            inputs = ModelParameterDictionary(input_stream)
        
        # Read input/configuration parameters
        self.m_n = inputs.get('MANNINGS_N', ptype=float)
        #NG do a check to make sure Manning's n is an appropriate value
        # below are default values, but it seems like in most cases we will not 
        # use default values
        self.rainfall_mmhr = inputs.get('RAINFALL_RATE', ptype=float)
        self.rain_duration = inputs.get('RAIN_DURATION', ptype=float)
        
        print "Rainfall duration ", self.rain_duration
        
        self.h_init = 0.001        # initial thin layer of water (m)
        self.g = 9.8               # gravitational acceleration (m/s2)
        self.alpha = 0.2           # time-step factor (ND; from Bates et al., 2010)
        self.m_n_sq = self.m_n*self.m_n # manning's n squared
        self.rho = 1000 # densith of water, kg/m^3
        
        
        # Derived parameters
        self.rainfall_rate = (self.rainfall_mmhr/1000.)/3600.  # rainfall in m/s
        self.ten_thirds = 10./3.   # pre-calculate 10/3 for speed
        
        # Set up state variables
        #NG should q be created outside of this Component?
        self.h = grid.create_node_dvector() + self.h_init     # water depth (m)
        print "length of h is ",self.h.size 
        #NG why is water depth at each node and not at cells, which are only active
        self.h_link = grid.create_active_link_dvector() #water depth at links (m)
        #this is needed for shear stress calculation
        self.q = grid.create_active_link_dvector()       # unit discharge (m2/s)
        self.dhdt = grid.create_active_cell_dvector()    # rate of water-depth change
        self.tau = grid.create_active_link_dvector()     # shear stress rho*g*h*S
        #NG why is dhdt at active cells, but h at nodes?
        #Maybe because you need a boundary condition of water depth at the 
        #boundary locations?  But water depths only chance at active cells?
        
    def run_one_step(self, grid, z, delt=None, rainrate=None):
        # In this case a delt should be a storm
        # and the units should be in seconds
        if delt==None:
            delt = self.rain_duration
        if rainrate==None:
            rainrate = self.rainfall_rate
        interior_cells = grid.get_active_cell_node_ids() 
        elapsed_time = 0
        
        while elapsed_time < delt:
            
            # Calculate time-step size for this iteration (Bates et al., eq 14)
            dtmax = self.alpha*grid.dx/np.sqrt(self.g*np.amax(self.h))
        
            # Take the smaller of delt or calculated time-step
            dt = min(dtmax, delt)
        
            #NG  need a loop in here that will cover the total time delt if 
            #dtmax < delt
        
            # Calculate the effective flow depth at active links. Bates et al. 2010
            # recommend using the difference between the highest water-surface
            # and the highest bed elevation between each pair of cells.
            zmax = grid.active_link_max(z) #array of length of num nodes
            w = self.h+z   # water-surface height
            wmax = grid.active_link_max(w) #array of length of num nodes
            hflow = wmax - zmax #array of length of num nodes
        
            # Calculate water-surface slopes: across links, but water heights are 
            #defined at nodes.
            water_surface_slope = grid.calculate_gradients_at_active_links(w)
       
            # Calculate the unit discharges (Bates et al., eq 11)
            self.q = (self.q-self.g*hflow*dt*water_surface_slope)/ \
                (1.+self.g*hflow*dt*self.m_n_sq*abs(self.q)/(hflow**self.ten_thirds))
            # q is calculated at links
            # water_surface_slope is at links
            # hflow is at links
            # calculate shear stress at links?  this seems crazy to NG, but give it
            # a go.

        
            # Calculate water-flux divergence at nodes
            dqds = grid.calculate_flux_divergence_at_nodes(self.q)
            
            # Second time-step limiter (experimental): make sure you don't allow
            # water-depth to go negative
            if np.amin(self.dhdt) < 0.:
                shallowing_locations = np.where(dhdt<0.)
                time_to_drain = -self.h[shallowing_locations]/self.dhdt[shallowing_locations]
                dtmax2 = alpha*np.amin(time_to_drain)
                dt = np.min([dtmax, dtmax2])
            else:
                dt = dtmax
        
            # Calculate rate of change of water depth
            self.dhdt = rainrate-dqds
        
            # Update the water-depth field
            self.h[interior_cells] = self.h[interior_cells] + self.dhdt[interior_cells]*dt
            #nic this needs to be fixed NICNICNIC
            #self.h[outlet_node] = self.h[node_next_to_outlet]

            # water_surface_slope is at links           
            # A problem is that the links have a direction, and you need to actually
            # know this direction when calculating shear stress.  
            # h is at nodes
            self.h_link = grid.assign_upslope_vals_to_active_links(self.h)
            #below is shear stress at the link
            #self.tau = self.rho*self.g*water_surface_slope*self.h_link
                                 
            # Update current time and return it
            self.current_time += dt
            
            elapsed_time += dt
            #taumax=np.amax(self.tau)
            #print "tau max ", taumax
            #watermax=np.amax(self.h)
            #print "water max ",watermax
            
            return self.tau
        
            #return z, g, qs, dqsds, dzdt
        
        
    #def run_one_step_internal(self, delt):
    #    
    #    # Take the smaller of delt or built-in time-step size self.dt
    #    dt = min(self.dt, delt)
    #    
    #    # Calculate the gradients and sediment fluxes
    #    self.g = self.grid.calculate_gradients_at_active_links(self.z)
    #    self.qs = -self.kd*self.g
    #    
    #    # Calculate the net deposition/erosion rate at each node
    #    self.dqsds = self.grid.calculate_flux_divergence_at_nodes(self.qs)
    #    
    #    # Calculate the total rate of elevation change
    #    dzdt = self.uplift_rate - self.dqsds
    #    
    #    # Update the elevations
    #    self.z[self.interior_cells] += dzdt[self.interior_cells] * dt
    #                             
    #    # Update current time and return it
    #    self.current_time += dt
    #    return self.current_time
    #    
    #    
    #def run_until_explicit(self, grid, z, t):
    #    
    #    while self.current_time < t:
    #        remaining_time = t - self.current_time
    #        z, g, qs, dqsds, dzdt = self.run_one_step_explicit(mg, z, g, qs, dqsds, dzdt, remaining_time)
    #        
    #    return z, g, qs, dqsds, dzdt
    #    
    #    
    #def run_until_internal(self, t):
    #    
    #    while self.current_time < t:
    #        remaining_time = t - self.current_time
    #        self.run_one_step_internal(remaining_time)
    #                
    #    
    #def get_time_step(self):
    #    """
    #    Returns time-step size.
    #    """
    #    return self.dt
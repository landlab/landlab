#! /usr/env/python
# -*- coding: utf-8 -*-
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
from matplotlib import pyplot as plt



class OverlandFlow(object):
    
    def __init__(self, input_stream, grid=None, current_time=0.):
        
        ## ALL CONSTANTS AND VARIABLES SHOULD BE DEFINED HERE. 
        
        #create and initial grid if one doesn't already exist

        self.grid = grid 
        if self.grid==None:
            self.grid = create_and_initialize_grid(input_stream)
        self.current_time = current_time
        self.initialize(grid, input_stream)

    def initialize(self, grid, input_stream, intensity=None, stormduration=None):
        
        # Create a ModelParameterDictionary for the inputs
        if type(input_stream)==ModelParameterDictionary:
            inputs = input_stream
        else:
            inputs = ModelParameterDictionary(input_stream)
        
        # Read input/configuration parameters
        self.m_n = inputs.get('MANNINGS_N', ptype=float)
        
        #NG do a check to make sure Manning's n is an appropriate values
        if intensity == None:
            self.rainfall_mmhr = inputs.get('RAINFALL_RATE', ptype=float)
        else:
            self.rainfall_mmhr = intensity
        
        if stormduration == None:
            self.rain_duration = inputs.get('RAIN_DURATION', ptype=float)
        else:
            self.rain_duration = stormduration
          
        self.h_init = 0.001        # initial thin layer of water (m)
        self.g = 9.8               # gravitational acceleration (m/s2)
        self.alpha = 0.2           # time-step factor (ND; from Bates et al., 2010)
        self.m_n_sq = self.m_n*self.m_n # manning's n squared
        self.rho = 1000 # densith of water, kg/m^3
              
        
        
        # Derived parameters
        self.rainfall_rate = (self.rainfall_mmhr/1000.)/3600.  # rainfall in m/s
        self.ten_thirds = 10./3.   # pre-calculate 10/3 for speed
        
        # Set up state variables
        
        self.hstart = grid.zeros(centering='node') + self.h_init 
        self.h = grid.zeros(centering='node') + self.h_init     # water depth (m)
        #NG why is water depth at each node and not at cells, which are only active
        self.q = grid.zeros(centering='active_link') # unit discharge (m2/s)
        self.dhdt = grid.zeros(centering='active_cell') # rate of water-depth change
        #NG why is dhdt at active cells, but h at nodes?
        #Maybe because you need a boundary condition of water depth at the 
        #boundary locations?  But water depths only change at active cells?
        
    def generate_overland_flow_at_one_point(self, grid, z, study_node, total_t=None, rainrate=None, rainduration=None):
        
        '''
        This method calculates discharge, water depth and shear stress at one
        point in the grid, defined as "study_node" in the function arguments.
        The study node is defined in the driver file, where the node ID must also
        be found using the grid_coords_to_node_id() function. 
        It is important to note that the study node does NOT have to be the outlet
        node, which must be defined to correctly account for boundary conditions.
        As it stands, using the outlet node as the "study_node" may cause significant
        boundary errors. Working on it!
        
        This function runs for the total run time, defined in the arguments as
        total_t (s). This function will work regardless of if rainduration (s) is shorter
        or longer than total_t (s)
        
        Rainrate is the rainfall intensity in m/s. 
        Rain duration is the total storm time in seconds.
        '''
    
        g=self.g
        alpha = self.alpha
        m_n_sq = self.m_n_sq
        ten_thirds = self.ten_thirds
        rho = self.rho
        
        if total_t==None:
            total_t = self.rain_duration
        if rainrate==None:
            rainrate = self.rainfall_rate
        if rainduration==None:
            rainduration = total_t
            
        #interior_nodes are the nodes on which you will be calculating flow 
        interior_nodes = grid.get_active_cell_node_ids()
        
        self.h = self.hstart
        h = self.h
        q = self.q
        dhdt = self.dhdt
        
        elapsed_time = 0
        
         # Get a list of the interior cells
        interior_cells = grid.get_active_cell_node_ids()
    
        # To track discharge at the outlet through time, we create initially empty
        # lists for time and outlet discharge.
        
        ## outlet_node AKA "study_node" = (input)
        q_outlet = []
        t = []
        t_1 = []
        h_1 = []
        q_outlet.append(0.)
        t.append(0.)
        t_1.append(0.)
        h_1.append(0.)
        
        nbr_node = grid.find_node_in_direction_of_max_slope_d4(z, study_node)
        study_link = grid.get_active_link_connecting_node_pair(study_node, nbr_node)


        # Main loop
        while elapsed_time < total_t:
        
        
            # Calculate time-step size for this iteration (Bates et al., eq 14)
            dtmax = alpha*grid.dx/np.sqrt(g*np.amax(h))
            #print "dtmax", dtmax
            
            if elapsed_time+dtmax > total_t:
                dtmax = total_t - elapsed_time
        
            # Calculate the effective flow depth at active links. Bates et al. 2010
            # recommend using the difference between the highest water-surface
            # and the highest bed elevation between each pair of cells.
            zmax = grid.max_of_link_end_node_values(z)
            w = h+z   # water-surface height
            wmax = grid.max_of_link_end_node_values(w)
            hflow = wmax - zmax
        
            # Calculate water-surface slopes
            water_surface_slope = grid.calculate_gradients_at_active_links(w)
       
            # Calculate the unit discharges (Bates et al., eq 11)
            q = (q-g*hflow*dtmax*water_surface_slope)/ \
                (1.+g*hflow*dtmax*0.06*0.06*abs(q)/(hflow**ten_thirds))
            
            #print "q study link", q[outlet_link] 
        
            # Calculate water-flux divergence at nodes
            dqds = grid.calculate_flux_divergence_at_nodes(q)
        
            # Update rainfall rate
            if elapsed_time > rainduration:
                rainrate = 0.
        
            # Calculate rate of change of water depth
            dhdt = rainrate-dqds
        
            # Second time-step limiter (experimental): make sure you don't allow
            # water-depth to go negative
            if np.amin(dhdt) < 0.:
                shallowing_locations = np.where(dhdt<0.)
                time_to_drain = -h[shallowing_locations]/dhdt[shallowing_locations]
                dtmax2 = alpha*np.amin(time_to_drain)
                dt = np.min([dtmax, dtmax2])
            else:
                dt = dtmax
        
            # Update the water-depth field
            h[interior_cells] = h[interior_cells] + dhdt[interior_cells]*dt

            w_slp_studypoint,garbage=grid.calculate_max_gradient_across_node_d4(w,study_node)
            tau_temp=self.rho*self.g*w_slp_studypoint*self.h[study_node]
            t_1.append(tau_temp)

            # Update model run time
            elapsed_time += dt
            print "elapsed time", elapsed_time
        
            # Remember discharge and time
            t.append(elapsed_time)
            q_outlet.append(q[study_link])            
            h_1.append(h[study_node])


        
        plt.figure('Discharge at Study Node')
        plt.plot(t, q_outlet, 'r-')
        plt.legend(loc=1)
        plt.ylabel('Discharge, m^3/s')
        plt.xlabel('Time, s')
        
        plt.figure('Shear Stress at Study Node')
        plt.plot(t, t_1, 'r--')
        plt.ylabel('Shear Stress, Pa')
        plt.xlabel('Time, s')
        plt.legend(loc=1)
        
        plt.figure('Water Depth at Study Node')
        plt.plot(t, h_1, 'r--')
        plt.ylabel('Water Depth, m')
        plt.xlabel('Time, s')
        plt.legend(loc=1)

     #   
        plt.show()
        #pylab.show()
        

    #def generate_overland_flow_across_grid()    
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

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

## PC RUNS

### TAU FILES ##
node1_t_file = 'C:\Users\Jordan\Dropbox\AGU RUNS\BaseCase\Trial1_5year\Tau\Node_1.txt'
node2_t_file = 'C:\Users\Jordan\Dropbox\AGU RUNS\BaseCase\Trial1_5year\Tau\Node_2.txt'
node3_t_file = 'C:\Users\Jordan\Dropbox\AGU RUNS\BaseCase\Trial1_5year\Tau\Node_3.txt'
#
### DISCHARGE FILES ##
node1_q_file = 'C:\Users\Jordan\Dropbox\AGU RUNS\BaseCase\Trial1_5year\Discharge\Node_1.txt'
node2_q_file = 'C:\Users\Jordan\Dropbox\AGU RUNS\BaseCase\Trial1_5year\Discharge\Node_2.txt'
node3_q_file = 'C:\Users\Jordan\Dropbox\AGU RUNS\BaseCase\Trial1_5year\Discharge\Node_3.txt'
#
### DEPTH FILES ##
node1_h_file = 'C:\Users\Jordan\Dropbox\AGU RUNS\BaseCase\Trial1_5year\Depth\Node_1.txt'
node2_h_file = 'C:\Users\Jordan\Dropbox\AGU RUNS\BaseCase\Trial1_5year\Depth\Node_2.txt'
node3_h_file = 'C:\Users\Jordan\Dropbox\AGU RUNS\BaseCase\Trial1_5year\Depth\Node_3.txt'



## MAC RUNS

### TAU FILES ##
#node1_t_file = '/Users/Jordan/Dropbox/AGU RUNS/BaseCase/Trial1_5year/Tau/Node_1.txt'
#node2_t_file = '/Users/Jordan/Dropbox/AGU RUNS/BaseCase/Trial1_5year/Tau/Node_2.txt'
#node3_t_file = '/Users/Jordan/Dropbox/AGU RUNS/BaseCase/Trial1_5year/Tau/Node_3.txt'
#
#### DISCHARGE FILES ##
#node1_q_file = '/Users/Jordan/Dropbox/AGU RUNS/BaseCase/Trial1_5year/Discharge/Node_1.txt'
#node2_q_file = '/Users/Jordan/Dropbox/AGU RUNS/BaseCase/Trial1_5year/Discharge/Node_2.txt'
#node3_q_file = '/Users/Jordan/Dropbox/AGU RUNS/BaseCase/Trial1_5year/Discharge/Node_3.txt'
#
#### DEPTH FILES ##
#node1_h_file = '/Users/Jordan/Dropbox/AGU RUNS/BaseCase/Trial1_5year/Depth/Node_1.txt'
#node2_h_file = '/Users/Jordan/Dropbox/AGU RUNS/BaseCase/Trial1_5year/Depth/Node_2.txt'
#node3_h_file = '/Users/Jordan/Dropbox/AGU RUNS/BaseCase/Trial1_5year/Depth/Node_3.txt'



def writetofile(namefile, arr):
    np.savetxt(namefile, arr,fmt ='%15.10f')

class OverlandFlow(object):
    
    def __init__(self, input_stream, grid=None, current_time=0.):
        
        self.grid = grid
        #create and initial grid if one doesn't already exist
        if self.grid==None:
            self.grid = create_and_initialize_grid(input_stream)
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
        #NG do a check to make sure Manning's n is an appropriate values
        self.rainfall_mmhr = inputs.get('RAINFALL_RATE', ptype=float)
        self.rain_duration = inputs.get('RAIN_DURATION', ptype=float)
        
        #print "Rainfall duration ", self.rain_duration
        
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
        
    def run_one_step(self, grid, z, outlet_node, delt=None, rainrate=None, rainduration=None):
        # run_one_step routes flow across the landscape over the time period delt.
        # This should work even if the rainduration is shorter or longer than delt. 
        # The units of delt and rainduration should be in seconds.
        # Note that the outlet_node and node_next_to_outlet do not need to be 
        # the actual outlet_node, just the location of interest for plotting.
        
        g=self.g
        alpha = self.alpha
        m_n_sq = self.m_n_sq
        ten_thirds = self.ten_thirds
        rho = self.rho
        
        if delt==None:
            delt = self.rain_duration
        if rainrate==None:
            rainrate = self.rainfall_rate
        if rainduration==None:
            rainduration = delt
            
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
        
        nbr_node = grid.find_node_in_direction_of_max_slope_d4(z, outlet_node)
        study_link = grid.get_active_link_connecting_node_pair(outlet_node, nbr_node)

                                 
        study_point2 = grid.grid_coords_to_node_id(138, 210)
        nbr_node2 = grid.find_node_in_direction_of_max_slope_d4(z, study_point2)
        study_link2 = grid.get_active_link_connecting_node_pair(study_point2, nbr_node2)
        h_2 = []
        h_2.append(0.)
        q_2 = []
        q_2.append(0.)
        t_2 = []
        t_2.append(0.)


        study_point3 = grid.grid_coords_to_node_id(140, 150)
        nbr_node3 = grid.find_node_in_direction_of_max_slope_d4(z, study_point3)
        study_link3 = grid.get_active_link_connecting_node_pair(study_point3, nbr_node3)
        h_3 = []
        h_3.append(0.)
        q_3 = []
        q_3.append(0.)
        t_3 = []
        t_3.append(0.)
        
        print "rainrate ",rainrate
        

        # Main loop
        while elapsed_time < delt:
        
        
            # Calculate time-step size for this iteration (Bates et al., eq 14)
            dtmax = alpha*grid.dx/np.sqrt(g*np.amax(h))
            #print "dtmax", dtmax
            
            if elapsed_time+dtmax > delt:
                dtmax = delt - elapsed_time
        
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

            w_slp_studypoint,garbage=grid.calculate_max_gradient_across_node_d4(w,outlet_node)
            tau_temp=self.rho*self.g*w_slp_studypoint*self.h[outlet_node]
            t_1.append(tau_temp)
                            
            w_slp_studypoint2,garbage=grid.calculate_max_gradient_across_node_d4(w,study_point2)
            tau_temp2=self.rho*self.g*w_slp_studypoint2*self.h[study_point2]
            t_2.append(tau_temp2)

            w_slp_studypoint3,garbage=grid.calculate_max_gradient_across_node_d4(w,study_point3)
            tau_temp3=self.rho*self.g*w_slp_studypoint3*self.h[study_point3]
            t_3.append(tau_temp3)

            # Update model run time
            elapsed_time += dt
            print "elapsed time", elapsed_time
        
            # Remember discharge and time
            t.append(elapsed_time)
            q_outlet.append(q[study_link])
            q_2.append(q[study_link2])
            q_3.append(q[study_link3])
            
            h_1.append(h[outlet_node])
            h_2.append(h[study_point2])
            h_3.append(h[study_point3])

        
        plt.figure('Discharge at Three Study Nodes')
        plt.plot(t, q_outlet, 'r-', label = 'NE Tributary')
        plt.plot(t, q_2, 'm-', label = 'E Tributary')
        plt.plot(t, q_3, 'c-', label = 'Upper Main Channel')
        plt.legend(loc=1)
        plt.ylabel('Discharge, m^3/s')
        plt.xlabel('Time, s')
        
        plt.figure('Shear Stress at Three Study Nodes')
        plt.plot(t, t_1, 'r--', label = 'NE Tributary')
        plt.plot(t, t_2, 'm--', label = 'E Tributary')
        plt.plot(t, t_3, 'c--', label = 'Upper Main Channel')
        plt.ylabel('Shear Stress, Pa')
        plt.xlabel('Time, s')
        plt.legend(loc=1)
        
        plt.figure('Water Depth at Three Study Nodes')
        plt.plot(t, h_1, 'r--', label = 'NE Tributary')
        plt.plot(t, h_2, 'm--', label = 'E Tributary')
        plt.plot(t, h_3, 'c--', label = 'Upper Main Channel')
        plt.ylabel('Water Depth, m')
        plt.xlabel('Time, s')
        plt.legend(loc=1)
        
        
        pylab.figure('Depth Map')
        hr = grid.node_vector_to_raster(h)
        im2 = pylab.imshow(hr, cmap=pylab.cm.PuBu,
                       extent=[0, grid.number_of_node_columns * grid.dx,
                               0, grid.number_of_node_rows * grid.dx])
        pylab.clim(0, 0.25)
        cb = pylab.colorbar(im2)
        cb.set_label('Water depth (m)', fontsize=12)
        pylab.title('Water depths for a 5 year storm')
        
        
        #plt.figure('Discharge at the NE Tributary')
        #plt.legend(loc=2)
        #plt.ylabel('Discharge, m^3/s')
        #plt.xlabel('Time, s')
        #plt.plot(t, q_outlet, 'r-', label = 'NE Tributary')
        #
        #plt.figure('Discharge at the E Tributary')
        #plt.plot(t, q_2, 'm-', label = 'E Tributary')
        #plt.legend(loc=2)
        #plt.ylabel("$Discharge, m^3/s$")
        #plt.xlabel('Time, s')
        #
     #   plt.figure('Discharge at the Upper Main Channel')
     #   plt.plot(t, q_3, 'c-', label = 'Upper Main Channel')
     #   plt.legend(loc=2)
     #   plt.ylabel('Discharge, m^3/s')
     #   plt.xlabel('Time, s')
     #
     #   plt.figure('Shear Stress at NE Tributary')
     #   plt.plot(t, t_1, 'r--', label = 'NE Tributary')
     #   plt.ylabel('Shear Stress, Pa')
     #   plt.xlabel('Time, s')
     #   plt.legend(loc=2)
     #   
     #   plt.figure('Shear Stress at E Tributary')
     #   plt.plot(t, t_2, 'm--', label = 'E Tributary')
     #   plt.ylabel('Shear Stress, Pa')
     #   plt.xlabel('Time, s')
     #   plt.legend(loc=2)
     #   
     #   plt.figure('Shear Stress at Upper Main Channel')
     #   plt.plot(t, t_3, 'c--', label = 'Upper Main Channel')
     #   plt.ylabel('Shear Stress, Pa')
     #   plt.xlabel('Time, s')
     #   plt.legend(loc=2)
     #   
     #   plt.figure('Water Depth at NE Tributary')
     #   plt.plot(t, h_1, 'r--', label = 'NE Tributary')
     #   plt.ylabel('Water Depth, m')
     #   plt.xlabel('Time, s')
     #   plt.legend(loc=2)     
     #   
     #   plt.figure('Water Depth at E Tributary')
     #   plt.plot(t, h_2, 'm--', label = 'E Tributary')
     #   plt.ylabel('Water Depth, m')
     #   plt.xlabel('Time, s')
     #   plt.legend(loc=2)
     #     
     #   plt.figure('Water Depth at Upper Main Channel')
     #   plt.plot(t, h_3, 'c--', label = 'Upper Main Channel')
     #   plt.ylabel('Water Depth, m')
     #   plt.xlabel('Time, s')
     #   plt.legend(loc=2)
        
        plt.show()
        pylab.show()
        
        writetofile(node1_t_file, t_1)
        writetofile(node2_t_file, t_2)
        writetofile(node3_t_file, t_3)

        writetofile(node1_h_file, h_1)
        writetofile(node2_h_file, h_2)
        writetofile(node3_h_file, h_3)

        writetofile(node1_q_file, q_outlet)
        writetofile(node2_q_file, q_2)
        writetofile(node3_q_file, q_3)    
        
        
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

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
        self.h = grid.zeros(centering='node') + self.h_init     # water depth (m)
        #NG why is water depth at each node and not at cells, which are only active
        #self.h_link = grid.zeros(centering='active_link') #water depth at links (m)
        #this is needed for shear stress calculation
        self.q = grid.zeros(centering='active_link') # unit discharge (m2/s)
        self.dhdt = grid.zeros(centering='active_cell') # rate of water-depth change
        #self.tau = grid.zeros(centering='active_link') # shear stress rho*g*h*S
        #NG why is dhdt at active cells, but h at nodes?
        #Maybe because you need a boundary condition of water depth at the 
        #boundary locations?  But water depths only change at active cells?
        
    def run_one_step(self, grid, z, outlet_node, node_next_to_outlet, delt=None, rainrate=None, rainduration=None):
        # In this case a delt is likely longer than a storm, in order to realize
        # the entire hydrograph.
        # The units should be in seconds
        if delt==None:
            delt = self.rain_duration
        if rainrate==None:
            rainrate = self.rainfall_rate
        if rainduration==None:
            rainduration = delt
        interior_nodes = grid.get_active_cell_node_ids()
        #print "interior nodes ", interior_nodes 
        #pylab.figure(5)
        #pylab.plot(np.array(interior_nodes),np.array(z[interior_nodes]))
        #pylab.ylabel('elevation of interior nodes')
        
        elapsed_time = 0
        
        #below is for calculating water surface slope at interior nodes
        w_slope=np.zeros(interior_nodes.size)
        #below is for calculating shear stress at interior nodes
        tau=np.zeros(interior_nodes.size)
        
        #print "length of w_slope ", w_slope.size   
        #print "delt ", delt
        
        
        #NG Everything below is for plotting purposes
        tau_plotter1 = grid.zeros(centering='node')
        h_plotter1 = grid.zeros(centering='node')
        
        plothelper=0
        #threetimes=[5,10,1000000000]
        threetimes=[239,799,1056,1000000000]
        
        study_point = grid.grid_coords_to_node_id(33,22)
        loc_study_point_interior = np.where(interior_nodes == study_point)
        h_study = []
        h_dwn = []
        tau_study = []
        q_study = []
        slope_study = []
        t = []
        h_dwn.append(0.)
        h_study.append(0.)
        q_study.append(0.)
        t.append(0.)
        tau_study.append(0.)
        slope_study.append(0.)
        study_link = grid.get_active_link_connecting_node_pair(outlet_node, 
                                                          node_next_to_outlet)
        helper=0
        #NG Done with plotting stuff.
        
        while elapsed_time < delt:
            
            # Calculate time-step size for this iteration (Bates et al., eq 14)
            dtmax = self.alpha*grid.dx/np.sqrt(self.g*np.amax(self.h))
            #print "dtmax", dtmax
        
            # Take the smaller of delt or calculated time-step
            dt = min(dtmax, delt)
            #print "dt", dt
        
            # Calculate the effective flow depth at active links. Bates et al. 2010
            # recommend using the difference between the highest water-surface
            # and the highest bed elevation between each pair of cells.
            zmax = grid.max_of_link_end_node_values(z) #array of length of active links
            w = self.h+z   # water-surface height, array of length num nodes
            wmax = grid.max_of_link_end_node_values(w) #array of length of active links
            hflow = wmax - zmax #array of length of activelinks
        
            # Calculate water-surface slopes: across links, but water heights are 
            #defined at nodes.
            water_surface_slope = grid.calculate_gradients_at_active_links(w)
       
            # Calculate the unit discharges (Bates et al., eq 11)
            self.q = (self.q-self.g*hflow*dt*water_surface_slope)/ \
                (1.+self.g*hflow*dt*self.m_n_sq*abs(self.q)/(hflow**self.ten_thirds))
            # q is calculated at links
            # water_surface_slope is at links
            # hflow is at links, but w is at nodes

            # Calculate water-flux divergence at nodes
            dqds = grid.calculate_flux_divergence_at_nodes(self.q)
            #dqds_max = np.amax(dqds)
            #print "max of dqds ", dqds_max
            
            # Calculate rate of change of water depth
            self.dhdt = rainrate-dqds
            #dhdt_max = np.amax(self.dhdt)
            #print "max of dhdt ", dhdt_max
            
            # Second time-step limiter (experimental): make sure you don't allow
            # water-depth to go negative
            if np.amin(self.dhdt) < 0.:
                shallowing_locations = np.where(self.dhdt<0.)
                time_to_drain = -self.h[shallowing_locations]/self.dhdt[shallowing_locations]
                dtmax2 = self.alpha*np.amin(time_to_drain)
                #print "dtmax2", dtmax2
                dt = np.min([dtmax, dtmax2, delt])
            #else:
            #   dt = np.min([dtmax,delt])
            #NG commented this out, because this was already calculated
        

            #watermax=np.amax(self.h)
            #print "water max before ",watermax
            #print "dt ", dt
        
            # Update the water-depth field
            self.h[interior_nodes] = self.h[interior_nodes] + self.dhdt[interior_nodes]*dt
            self.h[outlet_node] = self.h[node_next_to_outlet]

            # Let's calculate shear stress at the nodes.  
            # First get water height at the nodes.
            # Then calculate the maximum gradient in water surface elevations.
            # Then you can calculate shear stress!        
            # h is at nodes
            w = self.h+z   # water-surface height, array of length num nodes
            #print "length of w ", w.size
            
            #nbr_node = grid.find_node_in_direction_of_max_slope_d4(w, study_point)
            #study_link = grid.get_active_link_connecting_node_pair(study_point, 
            #                                              nbr_node)
                                                          
            #print "study node ",study_point," nbr node ", nbr_node, "other study point ", interior_nodes[loc_study_point_interior]
            #print "study link type ", type(study_link)
            
            #Below if is for limitingshear stress calculations to only times when
            #q surpasses a threshold
            if self.q[study_link]*grid.dx> 0.2:
                print "calculating shear stress, q study link ",self.q[study_link]*grid.dx
                for i in range(len(interior_nodes)): 
                    w_slope[i],garbage=grid.calculate_max_gradient_across_node_d4(w,interior_nodes[i])
                    tau[i]=self.rho*self.g*w_slope[i]*self.h[interior_nodes[i]]
                    #if interior_nodes[i] == study_point:
                    #    tau_study.append(tau[i])
                    #    slope_study.append(w_slope[i])
                        #print "found it"
                    
            #tau[np.where(tau<0)] = 0
                    
            #print "study slope ",w_slope[loc_study_point_interior]," study tau ",tau[loc_study_point_interior]
            #tau_study.append(tau[loc_study_point_interior])
            #slope_study.append(w_slope[loc_study_point_interior])        
            
            #OLD STUFF BELOW
            #self.h_link = grid.assign_upslope_vals_to_active_links(self.h)
            #below is shear stress at the link
            #self.tau = self.rho*self.g*water_surface_slope*self.h_link
                                 
            # Update current time and return it
            self.current_time += dt
            #print "elapsed_time ", elapsed_time
            elapsed_time += dt
            print "elapsed_time ", elapsed_time
            #taumax=np.max(tau)
            #taumin=np.min(tau)
            #print "tau max ", np.max(tau), " tau min ", np.min(tau)
            #print "water depth max ", np.max(self.h[interior_nodes]), " water depth min ", np.min(self.h[interior_nodes])
            #print "water surface slope ", np.max(w_slope), " water surface slope min ", np.min(w_slope)
            #print "water surface height ", np.max(w), " water surface height min ", np.min(w)
            if elapsed_time > rainduration:
                rainrate=0
                
            # Remember discharge and time
            t.append(elapsed_time)
            h_study.append(self.h[study_point])
            #h_dwn.append(self.h[nbr_node])
            
            #q is at links
            q_study.append(self.q[study_link])
            #print "q_study at ", helper, " is ", q_study[helper]
            #print "tau at ", helper, "is ", tau_study[helper]
            helper +=1
            
            #if elapsed_time > threetimes[plothelper]:
            #    for i in range(len(interior_nodes)): 
            #        tau_plotter1[interior_nodes[i]]=tau[i]
            #        h_plotter1[interior_nodes[i]]=self.h[i]
            #    plothelper +=1
            #    tr = grid.node_vector_to_raster(tau_plotter1)
            #    hr = grid.node_vector_to_raster(h_plotter1)
            #    pylab.figure(100)
            #    pylab.subplot(121)
            #    im = pylab.imshow(hr, cmap=pylab.cm.RdBu, extent=[0, grid.ncols*grid.dx, 0, grid.nrows*grid.dx])
            #    cb = pylab.colorbar(im)
            #    cb.set_label('flow depth (m)', fontsize=12)
            #    pylab.title('Flow Depth')
            #    pylab.subplot(122)
            #    im = pylab.imshow(tr, cmap=pylab.cm.RdBu, extent=[0, grid.ncols*grid.dx, 0, grid.nrows*grid.dx])
            #    cb = pylab.colorbar(im)
            #    cb.set_label('shear stress (Pa)', fontsize=12)
            #    pylab.title('Shear stress')
            #    pylab.show()
                
        #pylab.figure(1)
        ##pylab.plot(np.array(t), np.array(h_study))
        #pylab.plot(t,h_study)
        #pylab.xlabel('Time (s)')
        #pylab.ylabel('h (m)')
        #pylab.title('study point water height')
        #np.savetxt('h3322.data',h_study,fmt='%15.10f',delimiter='\n')
        #
        #np.savetxt('time.data',t,fmt='%15.10f',delimiter='\n')
        
        #pylab.figure(2)
        #pylab.plot(np.array(t), np.array(tau_study))
        #pylab.xlabel('Time (s)')
        #pylab.ylabel('shear stress (kg/m/s2)')
        #pylab.title('study point shear stress')
        #np.savetxt('tau3322.data',tau_study,fmt='%15.10f',delimiter='\n')
        
        # Plot discharge vs. time
        pylab.figure(3)
        pylab.plot(np.array(t), np.array(q_study)*grid.dx)
        #pylab.plot(t, q_study)
        pylab.xlabel('Time (s)')
        pylab.ylabel('Q (m3/s)')
        pylab.title('study point discharge')
        #
        #np.savetxt('q3322.data',q_study,fmt='%15.10f',delimiter='\n')
        
        #pylab.figure(4)
        #pylab.plot(np.array(t), np.array(slope_study))
        #pylab.xlabel('Time (s)')
        #pylab.ylabel('Water Slope (.)')
        #pylab.title('Water Slope at study point')
        #np.savetxt('slope3322.data',slope_study,fmt='%15.10f',delimiter='\n')
        
        #pylab.figure(5)
        #pylab.plot(np.array(t), np.array(h_dwn))
        #pylab.xlabel('Time (s)')
        #pylab.ylabel('h (m)')
        #pylab.title('down stream water height')
        #np.savetxt('hdwn3322.data',h_dwn,fmt='%15.10f',delimiter='\n')
        #
        # Display the plots
        pylab.show()
            
        
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

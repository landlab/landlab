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
        self.q = grid.zeros(centering='active_link') # unit discharge (m2/s)
        self.dhdt = grid.zeros(centering='active_cell') # rate of water-depth change
        #NG why is dhdt at active cells, but h at nodes?
        #Maybe because you need a boundary condition of water depth at the 
        #boundary locations?  But water depths only change at active cells?
        
    def run_one_step(self, grid, z, outlet_node, node_1, node_2, node_3, node_4, node_5, node_6, node_7, node_8, node_9, delt=None, rainrate=None, rainduration=None):
        # run_one_step routes flow across the landscape over the time period delt.
        # This should work even if the rainduration is shorter or longer than delt. 
        # The units of delt and rainduration should be in seconds.
        # Note that the outlet_node and node_next_to_outlet do not need to be 
        # the actual outlet_node, just the location of interest for plotting.
        
        if delt==None:
            delt = self.rain_duration
        if rainrate==None:
            rainrate = self.rainfall_rate
        if rainduration==None:
            rainduration = delt
            
        #interior_nodes are the nodes on which you will be calculating flow 
        #interior_nodes = grid.get_active_cell_node_ids()
        
        elapsed_time = 0
        
        #below is for calculating water surface slope at interior nodes
        #w_slope=np.zeros(interior_nodes.size)
        #below is for calculating shear stress at interior nodes
        #tau=np.zeros(interior_nodes.size)
        
        #print "length of w_slope ", w_slope.size   
        #print "delt ", delt
        
        
        #NG Everything below is for plotting purposes
        #tau_plotter1 = grid.zeros(centering='node')
        #above is for plotting shear stress at each nod
        #h_plotter1 = grid.zeros(centering='node')
        #above is for plotting flow depth at each node
        #
        #plothelper=0
        #threetimes=[239,799,1056,1000000000]
        #above is if you want to plot values at the first three times in the array
        
        study_point = outlet_node #redundant, but leave for now
        h_study = [] #flow depth at study node
        h_study.append(0.) #initialize array
        q_study = [] #discharge at study link
        q_study.append(0.) #initialize array
        #NG for shear stress, if calculating
        tau_study = [] #shear stress at study node
        tau_study.append(0.) #initialize array
        
        t = [] #time array for plotting
        t.append(0.) #initialize array
        
        study_point1 = grid.grid_coords_to_node_id(node_1)

        h_1 = []
        h_1.append(0.)
        q_1 = []
        q_1.append(0.)
        t_1 = []
        t_1.append(0.)
        
        study_point2 = grid.grid_coords_to_node_id(node_2)
 
        h_2 = []
        h_2.append(0.)
        q_2 = []
        q_2.append(0.)
        t_2 = []
        t_2.append(0.)
        
        study_point3 = grid.grid_coords_to_node_id(node_3)

        h_3 = []
        h_3.append(0.)
        q_3 = []
        q_3.append(0.)
        t_3 = []
        t_3.append(0.)
        
        study_point4 = grid.grid_coords_to_node_id(node_4)
        
        h_4 = []
        h_4.append(0.)
        q_4 = []
        q_4.append(0.)
        t_4 = []
        t_4.append(0.)
        
        study_point5 = grid.grid_coords_to_node_id(node_5)
        
        h_5 = []
        h_5.append(0.)
        q_5 = []
        q_5.append(0.)
        t_5 = []
        t_5.append(0.)
        
        study_point6 = grid.grid_coords_to_node_id(node_6)
        
        h_6 = []
        h_6.append(0.)
        q_6 = []
        q_6.append(0.)
        t_6 = []
        t_6.append(0.)
        
        study_point7 = grid.grid_coords_to_node_id(node_7)
        
        h_7 = []
        h_7.append(0.)
        q_7 = []
        q_7.append(0.)
        t_7 = []
        t_7.append(0.)
        
        study_point8 = grid.grid_coords_to_node_id(node_8)
        
        h_8 = []
        h_8.append(0.)
        q_8 = []
        q_8.append(0.)
        t_8 = []
        t_8.append(0.)
        
        study_point9 = grid.grid_coords_to_node_id(node_9)
        
        h_9 = []
        h_9.append(0.)
        q_9 = []
        q_9.append(0.)
        t_9 = []
        t_9.append(0.)
        
        #JORDAN, if you want to track values at another point, just edit below
        #to the correct coords.  Takes row then column.
        #study_point2 = grid.grid_coords_to_node_id(234,125)
        #h_study2 = [] #flow depth at study node 2
        #h_study2.append(0.) #initialize array 
        #q_study2 = [] #discharge at study link 2
        #q_study2.append(0.) #initialize array
        #tau_study2 = [] #shear stress at study node 2
        #tau_study2.append(0.) #initialize array    
        
        #NG for plotting purposes
        #slope_study = [] #water surface slope at study node
        #slope_study.append(0.) #initialize array
        
        #discharge is calculated at links, so you need the study link for finding
        #the discharge
        #study_link = grid.get_active_link_connecting_node_pair(outlet_node, 
        #                                                  node_next_to_outlet)
        #print "study link ",study_link
        #helper=0
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
            hflow = wmax - zmax #array of length of active links
        
            # Calculate water-surface slopes: across links, but water heights are 
            #defined at nodes.
            water_surface_slope = grid.calculate_gradients_at_active_links(w)
       
            # Calculate the unit discharges (Bates et al., eq 11)
            self.q = (self.q-self.g*hflow*dt*water_surface_slope)/ \
                (1.+self.g*hflow*dt*self.m_n_sq*abs(self.q)/(hflow**self.ten_thirds))
            
            #print "q study link", self.q[study_link]
                
            #NOTES:    
            # q is calculated at links
            # water_surface_slope is at links
            # hflow is at links, but w is at nodes

            # Calculate water-flux divergence at nodes
            dqds = grid.calculate_flux_divergence_at_nodes(self.q)
            
            # Calculate rate of change of water depth
            self.dhdt = rainrate-dqds
            
            #DEBUGGING
            #dqds_max = np.amax(dqds)
            #print "max of dqds ", dqds_max
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
            #NG commented out else, because this was already calculated

        
            # Update the water-depth field
            #self.h[interior_nodes] = self.h[interior_nodes] + self.dhdt[interior_nodes]*dt
            #self.h[outlet_node] = self.h[node_next_to_outlet]
            
            # Let's calculate shear stress at the nodes.  
            # First get water height at the nodes.
            # Then calculate the maximum gradient in water surface elevations (S).
            # Then you can calculate shear stress! (rho g h S)       
            # h (water depth) is at nodes
            
            w = self.h+z   # water-surface height, array of length num nodes
            
            #Below is a different way for finding the study discharge.
            #Rather than always plotting at a prescribed link, below will  
            #find the node next to the study point
            #that is connected by the steepest water surface slope.
            #This link is being used to calculate shear stress, so seems like 
            #a better discharge to use for plotting.
            #Note that this link was already set above, and so if this is 
            #commented out, code will still work.  The link above is likely
            #not the same as this link.
            
            nbr_node = grid.find_node_in_direction_of_max_slope_d4(w, study_point)
            study_link = grid.get_active_link_connecting_node_pair(study_point, nbr_node)
            
            nbr_node1 = grid.find_node_in_direction_of_max_slope_d4(w, study_point1)
            study_link1 = grid.get_active_link_connecting_node_pair(study_point1, nbr_node1)   
                        
            nbr_node2 = grid.find_node_in_direction_of_max_slope_d4(w, study_point2)
            study_link2 = grid.get_active_link_connecting_node_pair(study_point2, nbr_node2)      
            
            nbr_node3 = grid.find_node_in_direction_of_max_slope_d4(w, study_point3)
            study_link3 = grid.get_active_link_connecting_node_pair(study_point3, nbr_node3)   
            
            nbr_node4 = grid.find_node_in_direction_of_max_slope_d4(w, study_point4)
            study_link4 = grid.get_active_link_connecting_node_pair(study_point4, nbr_node4)      
            
            nbr_node5 = grid.find_node_in_direction_of_max_slope_d4(w, study_point5)
            study_link5 = grid.get_active_link_connecting_node_pair(study_point5, nbr_node5)   
            
            nbr_node6 = grid.find_node_in_direction_of_max_slope_d4(w, study_point6)
            study_link6 = grid.get_active_link_connecting_node_pair(study_point6, nbr_node6)      
            
            nbr_node7 = grid.find_node_in_direction_of_max_slope_d4(w, study_point7)
            study_link7 = grid.get_active_link_connecting_node_pair(study_point7, nbr_node7)   
            
            nbr_node8 = grid.find_node_in_direction_of_max_slope_d4(w, study_point8)
            study_link8 = grid.get_active_link_connecting_node_pair(study_point8, nbr_node8)      
            
            nbr_node9 = grid.find_node_in_direction_of_max_slope_d4(w, study_point9)
            study_link9 = grid.get_active_link_connecting_node_pair(study_point9, nbr_node9)   
            
            
                                    
                                                                    
            #print "study node ",study_point," nbr node ", nbr_node, "study_link ", study_link
            
            #JORDAN - uncomment below for shear stress calculations everywhere.
            #Below if is for limiting shear stress calculations to only times when
            #q surpasses a threshold (in this case q should be in m^3/sec)
            #Note that this threshold is hardwired below (on right of >)
            #This should be changed eventually.
            
            #if self.q[study_link]*grid.dx> 0.2:
            #    print "calculating shear stress, q study link ",self.q[study_link]*grid.dx
            #    for i in range(len(interior_nodes)): 
            #        w_slope[i],garbage=grid.calculate_max_gradient_across_node_d4(w,interior_nodes[i])
            #        tau[i]=self.rho*self.g*w_slope[i]*self.h[interior_nodes[i]]
            #        #if interior_nodes[i] == study_point:
            #        #    tau_study.append(tau[i])
            #        #    slope_study.append(w_slope[i])
            #            #print "found it"
                    
            #tau[np.where(tau<0)] = 0
                    
            #print "study slope ",w_slope[study_point]," study tau ",tau[study_point]
            #tau_study.append(tau[study_point])
            #slope_study.append(w_slope[study_point]) 
            
           #JORDAN, second attempt at calculating shear stress, but now I will 
           #calculate just at the study point, rather than at all points.
           
            w_slp_studypoint,garbage=grid.calculate_max_gradient_across_node_d4(w,study_point)
            tau_temp=self.rho*self.g*w_slp_studypoint*self.h[study_point]
            tau_study.append(tau_temp)

            w_slp_studypoint1,garbage=grid.calculate_max_gradient_across_node_d4(w,study_point1)
            tau_temp1=self.rho*self.g*w_slp_studypoint1*self.h[study_point1]
            t_1.append(tau_temp1)
                                                                  
            w_slp_studypoint2,garbage=grid.calculate_max_gradient_across_node_d4(w,study_point2)
            tau_temp2=self.rho*self.g*w_slp_studypoint2*self.h[study_point2]
            t_2.append(tau_temp2)

            w_slp_studypoint3,garbage=grid.calculate_max_gradient_across_node_d4(w,study_point3)
            tau_temp3=self.rho*self.g*w_slp_studypoint3*self.h[study_point3]
            t_3.append(tau_temp3)
                                                                  
            w_slp_studypoint4,garbage=grid.calculate_max_gradient_across_node_d4(w,study_point4)
            tau_temp4=self.rho*self.g*w_slp_studypoint4*self.h[study_point4]
            t_4.append(tau_temp4)
                                 
            w_slp_studypoint5,garbage=grid.calculate_max_gradient_across_node_d4(w,study_point5)
            tau_temp5=self.rho*self.g*w_slp_studypoint5*self.h[study_point5]
            t_5.append(tau_temp5)

            w_slp_studypoint6,garbage=grid.calculate_max_gradient_across_node_d4(w,study_point6)
            tau_temp6=self.rho*self.g*w_slp_studypoint6*self.h[study_point6]
            t_6.append(tau_temp6)

            w_slp_studypoint7,garbage=grid.calculate_max_gradient_across_node_d4(w,study_point7)
            tau_temp7=self.rho*self.g*w_slp_studypoint7*self.h[study_point7]
            t_7.append(tau_temp7)
                                 
            w_slp_studypoint8,garbage=grid.calculate_max_gradient_across_node_d4(w,study_point8)
            tau_temp8=self.rho*self.g*w_slp_studypoint8*self.h[study_point8]
            t_8.append(tau_temp8)

            w_slp_studypoint9,garbage=grid.calculate_max_gradient_across_node_d4(w,study_point9)
            tau_temp9=self.rho*self.g*w_slp_studypoint9*self.h[study_point9]
            t_9.append(tau_temp9)
                                                                                                                                                                                                                                                                                                                     
            # Update current time and return it
            #NG not sure what current_time is used for
            self.current_time += dt
            elapsed_time += dt
            print "elapsed_time ", elapsed_time
            
            if elapsed_time > rainduration:
                rainrate=0
                
            #NG used below for debugging    
            #taumax=np.max(tau)
            #taumin=np.min(tau)
            #print "tau max ", np.max(tau), " tau min ", np.min(tau)
            #print "water depth max ", np.max(self.h[interior_nodes]), " water depth min ", np.min(self.h[interior_nodes])
            #print "water surface slope ", np.max(w_slope), " water surface slope min ", np.min(w_slope)
            #print "water surface height ", np.max(w), " water surface height min ", np.min(w)
                
            # Remember discharge and time
            #Jordan, if you decide to track discharge and water depth at more 
            #than one point, you need to add some variable setting here.
            t.append(elapsed_time)
            h_study.append(self.h[study_point])
            h_1.append(self.h[study_point1])
            h_2.append(self.h[study_point2])
            h_3.append(self.h[study_point3])
            h_4.append(self.h[study_point4])
            h_5.append(self.h[study_point5])
            h_6.append(self.h[study_point6])
            h_7.append(self.h[study_point7])
            h_8.append(self.h[study_point8])
            h_9.append(self.h[study_point9])
            #h_dwn.append(self.h[nbr_node])
            
            #q is at links
            q_study.append(self.q[study_link])
            q_1.append(self.q[study_link1])
            q_2.append(self.q[study_link2])
            q_3.append(self.q[study_link3])
            q_4.append(self.q[study_link4])
            q_5.append(self.q[study_link5])
            q_6.append(self.q[study_link6])
            q_7.append(self.q[study_link7])
            q_8.append(self.q[study_link8])
            q_9.append(self.q[study_link9])
            
            
            
            #print "q_study is ", self.q[study_link]
            #print "tau at ", helper, "is ", tau_study[helper]
            #helper +=1
            
            #Below NG was using to plot data at three different times during
            #hydrograph
            #if elapsed_time > threetimes[plothelper]:
            #    for i in range(len(interior_nodes)): 
            #        tau_plotter1[interior_nodes[i]]=tau[i]
            #        h_plotter1[interior_nodes[i]]=self.h[i]
            #    plothelper +=1
            #    tr = grid.node_vector_to_raster(tau_plotter1)
            #    hr = grid.node_vector_to_raster(h_plotter1)
            #    pylab.figure(100)
            #    pylab.subplot(121)
            #    im = pylab.imshow(hr, cmap=pylab.cm.RdBu, extent=[0, grid.number_of_node_columns*grid.dx, 0, grid.number_of_node_rows*grid.dx])
            #    cb = pylab.colorbar(im)
            #    cb.set_label('flow depth (m)', fontsize=12)
            #    pylab.title('Flow Depth')
            #    pylab.subplot(122)
            #    im = pylab.imshow(tr, cmap=pylab.cm.RdBu, extent=[0, grid.number_of_node_columns*grid.dx, 0, grid.numer_of_node_rows*grid.dx])
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
        
        # Plot shear stress vs. time
        pylab.figure(4)
        pylab.plot(np.array(t), np.array(tau_study))
        pylab.xlabel('Time (s)')
        pylab.ylabel('shear stress (rho g h S)')
        pylab.title('study point shear stress')
        
        #Below for saving data to a text file
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

# -*- coding: utf-8 -*-
""" generate_overland_flow.py 

 This component simulates overland flow using
 the 2-D numerical model of shallow-water flow
 over topography using the Bates et al. (2010)
 algorithm for storage-cell inundation modeling.

 Most of code written by Greg Tucker.
 Shear stress component added by Nicole Gasparini


"""

from landlab import create_and_initialize_grid
import pylab
from landlab import ModelParameterDictionary
import numpy as np
from matplotlib import pyplot as plt
from math import atan, degrees
from landlab.components.flow_routing.route_flow_dn import FlowRouter
import os


_DEFAULT_INPUT_FILE = os.path.join(os.path.dirname(__file__),
                                  'input_data.txt')

def writetofile(namefile, arr):
    np.savetxt(namefile, arr,fmt ='%15.10f')
    
def change_signs(inarr1, outarr1):
    for i in inarr1:
        if i < 0:
            ind=inarr1.index(i)
            inarr1[ind] = -1
    for j in outarr1:
        if j > 0:
            ind2=outarr1.index(j)
            outarr1[ind2] = -1
    return inarr1, outarr1
    



class OverlandFlow(object):
    
    def __init__(self):
        """Simulation of overland flow using the shallow water equations
        discussed in Bates et al., (2010)
        """
        ## we should definitely reiterate that an initial depth is set (albiet small) -
        ## should probably be moved to the parameter text file...
        
        self.grid = None
        self.input_file = None

        self.current_time = 0.0
        self.m_n = 0.0
        self.intensity = 0.0
        self.stormduration = 0.0
        self.m_n_sq = self.m_n*self.m_n # manning's n squared
       
        self.h_init = 0.001           ## initial depth
        self.g = 9.8               # gravitational acceleration (m/s2)
        self.alpha = 0.2           # time-step factor (ND; from Bates et al., 2010)
        self.rho = 1000 # density of water, kg/m^3  
        self.ten_thirds = 10./3.   # pre-calculate 10/3 for speed



    def initialize(self, grid=None, input_file=None, intensity=None, stormduration=None):
        """ Initialize the model and grid using either an input file
        or default values. 
        """
        self.grid = grid 
        if self.grid==None:
            self.grid = create_and_initialize_grid(input_file) ##<- this is the same input file used for parameters. 
                                                                ## this seems wrong to me...
            
        ##self.current_time = current_time <- this isn't being used yet. 
        
        # Create a ModelParameterDictionary for the inputs
        
        MPD = ModelParameterDictionary()
        ## I am not sure we want to use an instance of MPD as an input_file. I'm confused about this. SO for 
        ## now I changed this to reflect other component types
        
        if input_file is None:
            input_file = _DEFAULT_INPUT_FILE
            
        MPD.read_from_file(input_file)
        
        ## This was the old way... where an instance of MPD is an input_file
        #if type(input_file)==MPD: 
        #    inputs = input_file
        #else:
        #    inputs = MPD(input_file)
        
        # Read input/configuration parameters
        self.m_n = MPD.read_float('MANNINGS_N')
        
        if intensity == None:
            self.rainfall_mmhr = MPD.read_float( 'RAINFALL_RATE')
        else:
            self.rainfall_mmhr = intensity
        
        if stormduration == None:
            self.rain_duration = MPD.read_float( 'RAIN_DURATION' )
        else:
            self.rain_duration = stormduration
          
        #self.h_init = 0.001        # initial thin layer of water (m)
        #self.g = 9.8               # gravitational acceleration (m/s2)
        #self.alpha = 0.2           # time-step factor (ND; from Bates et al., 2010)
        #self.m_n_sq = self.m_n*self.m_n # manning's n squared
        #self.rho = 1000 # density of water, kg/m^3     
        
        # Derived parameters ## THIS IS IMPORTANT - UNITS!!!! Double check precip component...
        ## IF WE CLEARLY STATE UNIT REQUIREMENTS IN DOCS, THIS SHOULDN'T BE A PROBLEM
        ## BUT IT MUST BE CLEAR THAT ALL INPUT FILES HAVE STANDARD UNITS THROUGHOUT. 
        
        self.rainfall_rate = (self.rainfall_mmhr/1000.)/3600.  # rainfall in m/s
        
        
        # Set up state variables
        
        self.hstart = grid.zeros(centering='node') + self.h_init 
        self.h = grid.zeros(centering='node') + self.h_init     # water depth (m)
        #NG why is water depth at each node and not at cells, which are only active
        self.q = grid.zeros(centering='active_link') # unit discharge (m2/s)
        self.dhdt = grid.zeros(centering='cell') # rate of water-depth change
        
        self.tau = grid.zeros(centering='node') # shear stress (Pascals)
        #self.dtaudt = grid.zeros(centering='active_cell') # rate of shear stress change
        #NG why is dhdt at active cells, but h at nodes?
        #Maybe because you need a boundary condition of water depth at the 
        #boundary locations?  But water depths only change at active cells?
        
    def flow_at_one_node(self, grid, z, study_node, total_t=None, rainrate=None, rainduration=None):
        
        '''This method calculates discharge, water depth and shear stress at one
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
            

        self.h = self.hstart
        h = self.h
        q = self.q
        dhdt = self.dhdt
        
        elapsed_time = 0
        
         # Get a list of the interior cells
        self.interior_nodes = grid.get_active_cell_node_ids()

    
        # To track discharge at the study node through time, we create initially empty
        # lists for the node discharge, tau and depth. Time is also saved to help with
        # plotting later...
        
        ## outlet_node AKA "study_node" = (input)
        self.q_study = []
        self.time = []
        self.tau_study = []
        self.depth_study = []
        self.q_study.append(0.)
        self.time.append(0.)
        self.tau_study.append(0.)
        self.depth_study.append(0.)
        
        nbr_node = grid.find_node_in_direction_of_max_slope_d4(z, study_node)
        study_link = grid.get_active_link_connecting_node_pair(study_node, nbr_node)
        self.q_node = grid.zeros(centering='node')#grid.get_active_cell_node_ids()
        self.tau_node = grid.zeros(centering='node')
        self.dqds = grid.zeros(centering='node')
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
            self.w = h+z   # water-surface height
            wmax = grid.max_of_link_end_node_values(self.w)
            hflow = wmax - zmax
        
            # Calculate water-surface slopes
            water_surface_slope = grid.calculate_gradients_at_active_links(self.w)
       
            # Calculate the unit discharges (Bates et al., eq 11)
            self.q = (self.q-g*hflow*dtmax*water_surface_slope)/ \
                (1.+g*hflow*dtmax*0.06*0.06*abs(self.q)/(hflow**ten_thirds))
            
            #print "q study link", q[outlet_link] 
        
            # Calculate water-flux divergence at nodes
            self.dqds = grid.calculate_flux_divergence_at_nodes(self.q)
            #dqds = grid.calculate_gradients_at_active_links(q)
        
            # Update rainfall rate
            if elapsed_time > rainduration:
                rainrate = 0.
        
            # Calculate rate of change of water depth
            dhdt = rainrate-self.dqds
        
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
            h[self.interior_nodes] = h[self.interior_nodes] + dhdt[self.interior_nodes]*dt

            self.slopes_at_node, garbage = grid.calculate_steepest_descent_on_nodes(self.w, water_surface_slope)
            w_slp_studypoint,garbage=grid.calculate_max_gradient_across_node_d4(self.w,study_node)
            tau_temp=self.rho*self.g*w_slp_studypoint*self.h[study_node]
            self.tau_study.append(tau_temp)
            self.depth_study.append(h[study_node])

            # Update model run time
            elapsed_time += dt
            print "elapsed time", elapsed_time
        
            # Remember discharge and time
            self.time.append(elapsed_time)
            self.q_study.append(q[study_link])            
            self.depth_study.append(h[study_node])

    def plot_at_one_node(self):
        
        '''This method must follow a call to the flow_at_one_node() method.
        It plots depth, discharge and shear stress through time at the study
        node. 
        '''       
        plt.figure('Discharge at Study Node')
        plt.plot(self.time, self.q_study, 'r-')
        plt.legend(loc=1)
        plt.ylabel('Discharge, m^3/s')
        plt.xlabel('Time, s')
        
        plt.figure('Shear Stress at Study Node')
        plt.plot(self.time, self.tau_study, 'r--')
        plt.ylabel('Shear Stress, Pa')
        plt.xlabel('Time, s')
        plt.legend(loc=1)
        
        plt.figure('Water Depth at Study Node')
        plt.plot(self.time, self.depth_study, 'r--')
        plt.ylabel('Water Depth, m')
        plt.xlabel('Time, s')
        plt.legend(loc=1)
        plt.show()
        
    def plot_water_depths(self, grid):
         ''' This plots water depths across the raster. '''
         plt.figure('Water Depth Raster')
         hr = grid.node_vector_to_raster(self.h)
         palette = pylab.cm.winter_r
         palette.set_under('w', 0.01) ## right now these are hardcoded values. 
         im2 = pylab.imshow(hr, cmap=palette, extent=[0, grid.number_of_node_columns *grid.dx,0, grid.number_of_node_rows * grid.dx])
         pylab.clim(0.01, 1) ## hardcoded - need to think of something clever.
         cb = pylab.colorbar(im2)
         cb.set_label('Water depth (m)', fontsize=12)
         pylab.title('Water depth')
         plt.show()
         
    def plot_discharge(self,grid):
        """This plots discharge across the raster"""
     
        grid._setup_active_inlink_and_outlink_matrices()
        outlink = grid.node_active_outlink_matrix
        inlink = grid.node_active_inlink_matrix
        outlink1 = outlink.tolist()
        inlink1 =inlink.tolist()
        newin0, newout0 = change_signs(inlink1[0], outlink1[0])
        newin1, newout1 = change_signs(inlink1[1], outlink1[1])
        in0 = np.array(newin0)
        in1 = np.array(newin1)
        out0 = np.array(newout0)
        out1 = np.array(newout1)
        self.q_node = self.q[in0]+self.q[in1]+self.q[out0]+self.q[out1] #((q[outlink[1]] -q[inlink[1]]))#+((q[outlink[0]]-(q[inlink[0]])))
        fixed_q = grid.zeros(centering='node')
        for each in self.interior_nodes:
            fixed_q[each] = self.q_node[each]
        plt.figure('DISCHARGE')
        hr = grid.node_vector_to_raster(fixed_q)
        palette = pylab.cm.RdYlBu
        im2 = pylab.imshow(hr, cmap=palette, extent=[0, grid.number_of_node_columns *grid.dx,0, grid.number_of_node_rows * grid.dx])
        pylab.clim(vmin=0.000001)#, mx)
        palette.set_under('w', 0.000001)
        cb = pylab.colorbar(im2)
        cb.set_label('DISCHARGE (m)', fontsize=12)
        pylab.title('DISCHARGE')
        plt.show()

    def plot_shear_stress_grid(self, grid):                            
        plt.figure('Shear Stress raster')
        tn = grid.node_vector_to_raster(self.tau)
        cmap=plt.get_cmap('Spectral', 10)
        cmap.set_under('white')
        im2 = pylab.imshow(tn, cmap=cmap, extent=[0, grid.number_of_node_columns *grid.dx,0, grid.number_of_node_rows * grid.dx])
        cb = pylab.colorbar(im2)
        pylab.clim(vmin=10, vmax=100)
        cb.set_label('Shear Stress (Pa)', fontsize=12)
        pylab.title('Shear Stress')
        plt.show()    
        
    def plot_slopes(self, grid):
        
        plt.figure('Water Surface Slopes')
        cmap=plt.get_cmap('RdYlGn_r', 10)
        fixed_slope = grid.zeros(centering='node')
        for each in self.interior_nodes:
            fixed_slope[each] = degrees(atan(self.slopes_at_node[each]))
        fx = grid.node_vector_to_raster(fixed_slope)
        cmap.set_under('w', 0.01) ## right now these are hardcoded values. 
        im2 = pylab.imshow(fx, cmap=cmap, extent=[0, grid.number_of_node_columns *grid.dx,0, grid.number_of_node_rows * grid.dx])
        pylab.clim(0.01, 50) ## hardcoded - need to think of something clever.
        cb = pylab.colorbar(im2)
        cb.set_label('Slope', fontsize=12)
        pylab.title('Slope')
        plt.show()
                    

    def flow_across_grid(self, grid, z, total_t=None, rainrate=None, rainduration=None):  
        '''
        This function calculates depth, discharge and shear stress
        at each active node in the grid. 
    
        As of right now, this is going to be incredibly time consuming.
        WORK IN PROGRESS!!!!!
        '''
        
        self.h = self.hstart
        h = self.h
        q = self.q
        dhdt = self.dhdt
        
        #self.tau
#        dtdt = self.dtaudt
        
        g=self.g
        alpha = self.alpha
        m_n_sq = self.m_n_sq
        ten_thirds = self.ten_thirds
        rho = self.rho
    
        #interior_nodes are the nodes on which you will be calculating flow 
        self.interior_nodes = grid.get_active_cell_node_ids()

        if total_t==None:
            total_t = self.rain_duration
        if rainrate==None:
            rainrate = self.rainfall_rate
        if rainduration==None:
            rainduration = total_t
            
        elapsed_time = 0
        
        #below is for calculating water surface slope at interior nodes
        w_slope=np.zeros(self.interior_nodes.size)

        t = [] #time array for plotting
        t.append(0.) #initialize array
        self.dqds= grid.zeros(centering='node')
                
        while elapsed_time < total_t:
            
            # Calculate time-step size for this iteration (Bates et al., eq 14)
            dtmax = self.alpha*grid.dx/np.sqrt(self.g*np.amax(self.h))
        
            # Take the smaller of delt or calculated time-step
            dt = min(dtmax, total_t)
        
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
            self.q = (self.q-g*hflow*dtmax*water_surface_slope)/ \
                (1.+g*hflow*dtmax*0.06*0.06*abs(self.q)/(hflow**ten_thirds))                
            #NOTES:    
            # q is calculated at links
            # water_surface_slope is at links
            # hflow is at links, but w is at nodes

            # Calculate water-flux divergence at nodes
            self.dqds = grid.calculate_flux_divergence_at_nodes(self.q,self.dqds)
            
            # Calculate rate of change of water depth
            self.dhdt = rainrate-self.dqds

            
            # Second time-step limiter (experimental): make sure you don't allow
            # water-depth to go negative
            if np.amin(self.dhdt) < 0.:
                shallowing_locations = np.where(self.dhdt<0.)
                time_to_drain = -self.h[shallowing_locations]/self.dhdt[shallowing_locations]
                dtmax2 = self.alpha*np.amin(time_to_drain)
                dt = np.min([dtmax, dtmax2, total_t])
        
            # Update the water-depth field
            self.h[self.interior_nodes] = self.h[self.interior_nodes] + self.dhdt[self.interior_nodes]*dt


            # Let's calculate shear stress at the nodes.  
            # First get water height at the nodes.
            # Then calculate the maximum gradient in water surface elevations (S).
            # Then you can calculate shear stress! (rho g h S)       
            # h (water depth) is at nodes
            
            w = self.h+z   # water-surface height, array of length num nodes
            

            #Below if is for limiting shear stress calculations to only times when
            #q surpasses a threshold (in this case q should be in m^3/sec)
            #Note that this threshold is HARDWIRED below (on right of >)

            #self.slopes_at_node, garbage = grid.calculate_steepest_descent_on_nodes(w, water_surface_slope)
            fr=FlowRouter(grid)
            r, a, q, ss, s, d = fr.route_flow(w)
            self.slopes_at_node = ss
            if self.q.any()*grid.dx> 0.2:
                self.tau = self.rho*self.g*self.slopes_at_node*self.h

                #for i in range(24899): 
                #    #w_slope[i],garbage=grid.calculate_gradient_across_cell_faces(grid, w, self.interior_nodes)    #grid.calculate_max_gradient_across_node_d4(w,self.interior_nodes[i])
                #    self.tau[self.interior_nodes[i]]=self.rho*self.g*self.slopes_at_node[self.interior_nodes[i]]*self.h[self.interior_nodes[i]]
                #    #self.tau[i] = tau[i]

            self.tau[np.where(self.tau<0)] = 0       

            #self.tau[interior_nodes] = self.tau[interior_nodes] + self.dtds[interior_nodes]*dt
            #tau[np.where(tau<0)] = 0
            # Update model run time
            elapsed_time += dt
            print "elapsed time", elapsed_time
            #print max(self.tau)
            #print min(self.tau)
            #print len(interior_nodes)
            

            #print "study slope ",w_slope[study_point]," study tau ",tau[study_point]
            #tau_study.append(tau[study_point])
            #slope_study.append(w_slope[study_point]) 
        
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

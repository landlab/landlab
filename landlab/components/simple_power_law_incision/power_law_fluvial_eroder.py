#! /usr/env/python
# -*- coding: utf-8 -*-
"""
Component for detachment-limited fluvial incision using a simple power-law model.

I = K Q^m S^n

I=incision rate (M Y^(-1) )
K=bedrock erodibility (M^(1-3m) Y^(m-1) )
Q= fluvial discharge (M^3 Y^-1 )
S=slope of landscape (negative of the gradient in topography, dimensionless, 
and only applies on positive slopes)

NOTE, in units above M=meters.  This component assumes that variables are given
in units of meters and years, including rainfall!

NOTE, Only incision happens in this class.  NO DEPOSITION.  
NO TRACKING OF SEDIMENT FLUX.

This assumes that a grid has already been instantiated.
 
"""

#import landlab
from landlab import ModelParameterDictionary
from landlab.components.flow_routing.flow_routing_D8 import RouteFlowD8
from landlab.components.flow_accum.flow_accumulation2 import AccumFlow
import numpy as np
import matplotlib.pyplot as plt

class PowerLawIncision(object):
    
    def __init__(self, input_stream, grid, current_time=0.):
        
        self.grid = grid
        #create and initial grid if one doesn't already exist
        #if self.grid==None:
        #    self.grid = create_and_initialize_grid(input_stream)
        
        self.current_time = current_time
        self.initialize(grid, input_stream)
        
    def initialize(self, grid, input_stream):
        
        # Create a ModelParameterDictionary for the inputs
        if type(input_stream)==ModelParameterDictionary:
            inputs = input_stream
        else:
            inputs = ModelParameterDictionary(input_stream)
        
        # Read input/configuration parameters
        self.m = inputs.get('M_EXPONENT', ptype=float)
        self.n = inputs.get('N_EXPONENT', ptype=float)
        self.K = inputs.get('K_COEFFICIENT', ptype=float)
        self.rainfall_myr = inputs.get('RAINFALL_RATE_M_PER_YEAR', ptype=float)
        self.rain_duration_yr = inputs.get('RAIN_DURATION_YEARS', ptype=float)
        self.frac = 0.9 #for time step calculations
        
        #print "Rainfall duration ", self.rain_duration
        
        # Set up state variables
        self.I = grid.zeros(centering='node')    # incision rate (M/Y)
        
    def run_one_storm(self, grid, z, rainrate=None, storm_dur=None):
        
        if rainrate==None:
            rainrate = self.rainfall_myr
        if storm_dur==None:
            storm_dur = self.rain_duration_yr
            
        m=self.m
        n=self.n
        K=self.K
        frac = self.frac
        
        #interior_nodes are the nodes on which you will be calculating incision 
        interior_nodes = grid.get_active_cell_node_ids()
        
        #instantiate variable of type RouteFlowD8 Class
        flow_router = RouteFlowD8(len(z))
        #initial flow direction
        flowdirs, max_slopes = flow_router.calc_flowdirs(grid,z)
        #insantiate variable of type AccumFlow Class
        accumulator = AccumFlow(grid)
        #initial flow accumulation
        drain_area = accumulator.calc_flowacc(grid, z, flowdirs)
        
        time=0
        dt = storm_dur
        while time < storm_dur:
            #Calculate incision rate, should be in m/yr, should be negative
            #First make sure that there are no negative (uphill slopes)
            #Set those to zero, because incision rate should be zero there.
            max_slopes.clip(0)
            I=-K*np.power(rainrate*drain_area,m)*np.power(max_slopes,n)
            
            #Do a time-step check
            #If the downstream node is eroding at a slower rate than the 
            #upstream node, there is a possibility of flow direction reversal,
            #or at least a flattening of the landscape. 
            #Limit dt so that this flattening or reversal doesn't happen.
            #How close you allow these two points to get to eachother is
            #determined by the variable frac.
            for i in interior_nodes:
                dzdtdif = I[flowdirs[i]]-I[i]
                if dzdtdif > 0:
                    dtmin = frac*(z[i]-z[flowdirs[i]])/dzdtdif
                    #Nic do you need a test here to make sure that dtmin 
                    #isn't too small?  Trying that out
                    if dtmin < dt:
                        if dtmin>0.001*dt:
                            dt = dtmin
                        else:
                            dt = 0.001*dt
            
            #should now have a stable timestep.
            #reduce elevations.
            z=I*dt+z
            
            #update elapsed time
            time=dt+time
            
            #check to see that you are within 0.01% of time
            #if so, done
            #otherwise, reset everything for next loop
            
            if time > 0.9999*storm_dur:
                #done!
                time = storm_dur
            else:
                #not done, reset everything
                #update time step to maximum possible
                dt = storm_dur - time
                #recalculate flow directions
                flowdirs, max_slopes = flow_router.calc_flowdirs(grid,z)
                #recalculate drainage area
                drain_area = accumulator.calc_flowacc(grid, z, flowdirs)
        
        return z       
        
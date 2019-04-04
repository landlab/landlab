#! /usr/env/python
# -*- coding: utf-8 -*-
"""
April 4, 2019 Starting to rehab lateral erosion
ALangston


"""

#import landlab
from landlab import ModelParameterDictionary
#from landlab.components.flow_routing.flow_routing_D8 import RouteFlowD8
from landlab.components.flow_routing import FlowRouter
from landlab.components.lateral_erosion.node_finder2 import Node_Finder2
#from landlab.components.radius_curv_dz import radius_curv_dz
#from landlab.components.flow_accum.flow_accumulation2 import AccumFlow
from landlab.utils import structured_grid
import numpy as np
np.set_printoptions(threshold=np.nan)
from random import uniform
import matplotlib.pyplot as plt

class LateralEroder(component):

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
        self.alph = inputs.get('ALPH', ptype=float)
        self.Kv = inputs.get('KV_COEFFICIENT', ptype=float)
        self.Klr = inputs.get('KL_RATIO', ptype=float)
        self.rain_duration_yr = inputs.get('RAIN_DURATION_YEARS', ptype=float)
        self.inlet_node = inputs.get('INLET_NODE', ptype=float)
        self.inlet_area = inputs.get('INLET_AREA', ptype=float)
        self.qsinlet = inputs.get('QSINLET', ptype=float)
        self.frac = 0.3 #for time step calculations

        # Set up state variables

        #initialize qsin for each interior node, all zero initially.
        self.qsin = grid.zeros(centering='node')    # qsin (M^3/Y)
        self.qt = grid.zeros(centering='node')    # transport capacity
        self.dzdt = grid.zeros(centering='node')    # elevation change rate (M/Y)
        self.qsqt = grid.zeros(centering='node')    # potential elevation change rate (M/Y)
    def save_multipagepdf(f_handle):    
    	savefig(f_handle, format='pdf')
    	close()


    def run_one_storm(self, grid, z, vol_lat, rainrate=None, storm_dur=None, qsinlet=None, inlet_area=None, Klr=None):

        if rainrate==None:
            rainrate = self.rainfall_myr
        if storm_dur==None:
            storm_dur = self.rain_duration_yr   
        inlet_node=self.inlet_node
        #inlet_area=self.inlet_area
        if qsinlet==None:
            qsinlet=self.qsinlet
        if inlet_area==None:
            inlet_area=self.inlet_area
        if Klr==None:    #Added10/9 to allow changing rainrate (indirectly this way.)
            Klr=self.Klr


        Kv=self.Kv
        #Klr=self.Klr
        frac = self.frac
        qsin=self.qsin
        qt=self.qt
        dzdt=self.dzdt
        qsqt=self.qsqt
        alph=self.alph
        #**********ADDED FOR WATER DEPTH CHANGE***************
        #now KL/KV ratio is a parameter set from infile again.
        #still need to calculate runoff for Q and water depth calculation
        kw=10.
        F=0.02
        #May 2, runoff calculated below (in m/s) is important for calculating
        #discharge and water depth correctly. renamed runoffms to prevent
        #confusion with other uses of runoff
        runoffms=(Klr*F/kw)**2

        #Kl is calculated from ratio of lateral to vertical K parameters
        Kl=Kv*Klr

        dx=grid.dx
        nr=grid.number_of_node_rows
        nc=grid.number_of_node_columns
        interior_nodes = grid.core_nodes
        boundary_nodes=structured_grid.boundary_nodes((nr,nc))
        #clear qsin for next loop
        qsin = grid.zeros(centering='node')
        qsqt = grid.zeros(centering='node')
        #eronode=np.zeros(grid.number_of_nodes)
        lat_nodes=np.zeros(grid.number_of_nodes)
        dzlat=np.zeros(grid.number_of_nodes)
        dzver=np.zeros(grid.number_of_nodes)
        vol_lat_dt=np.zeros(grid.number_of_nodes)
        #z_bank=np.zeros(grid.number_of_nodes)
	#vol_diff=np.zeros(grid.number_of_nodes)
        #instantiate variable of type RouteFlowDNClass
        flow_router = FlowRouter(grid)

        # 4/24/2017 add inlet to change drainage area with spatially variable runoff rate
        #runoff is an array with values of the area of each node (dx**2)
        runoffinlet=np.ones(grid.number_of_nodes)*dx**2
        #Change the runoff at the inlet node to node area + inlet node
        runoffinlet[inlet_node]=+inlet_area
        _=grid.add_field('node', 'water__unit_flux_in', runoffinlet,
                     noclobber=False)
        
        
        
#        flowdirs, drain_area, q, max_slopes, s, receiver_link = flow_router.route_flow(elevs=z, node_cell_area=node_area, runoff_rate=runoff)
        flow_router.route_flow(method='D8')
        #flow__upstream_node_order is node array contianing downstream to upstream order list of node ids
        s=grid.at_node['flow__upstream_node_order']
        drain_area_fr=grid.at_node['drainage_area']    #renamed this drainage area set by flow router
        max_slopes=grid.at_node['topographic__steepest_slope']
        q=grid.at_node['surface_water__discharge']
        flowdirs=grid.at_node['flow__receiver_node']
        drain_area=q/dx**2    #this is the drainage area that I need for code below with an inlet set by spatially varible runoff. 
        if (0):
#            print 'nodeIDs', grid.core_nodes
#            print 'flowupstream order', s
            print 'q', q.reshape(nr,nc)
            print 'runoffms', runoffms
            print 'q_ms', (drain_area*runoffms).reshape(nr,nc)
#            print 'drainareafr', drain_area_fr.reshape(nr,nc)
            print 'drain_area', drain_area.reshape(nr,nc)
            print delta
#        max_slopes2=grid.calc_grad_at_active_link(z)
#        max_slopes=grid.calc_slope_at_node()
#        print 'lengthflowdirs', len(flowdirs)
        #temporary hack for drainage area
#        drain_area[inlet_node]=inlet_area
        
        #order interior nodes
        #find interior nodes in downstream ordered vector s

        #make a list l, where node status is interior (signified by label 0) in s
        l=s[np.where((grid.status_at_node[s] == 0))[0]]
#        print 'l', l
        #this misses an interior nodes that is set as constant value, 1
        #but the below grabs nodes that are set as open boundaries. no good for me here
        #l2=s[np.where((grid.node_status[s] == 1))[0]]
        #combine lists
        #dwnst_nodes=np.insert(l,0,l2)
        dwnst_nodes=l
        #reverse list so we go from upstream to down stream
        #4/20/2017: this works
#        print "dwnst_nodes before reversal", dwnst_nodes
        dwnst_nodes=dwnst_nodes[::-1]
#        print "dwnst_nodes after reversal", dwnst_nodes
        #local time
        time=0
        dt = storm_dur
        dtmax=storm_dur

              

        while time < storm_dur:
            #Calculate incision rate, should be in m/yr, should be negative
            #First make sure that there are no negative (uphill slopes)
            #Set those to zero, because incision rate should be zero there.
#            print 'slopes', max_slopes
            max_slopes=max_slopes.clip(0)
#            print 'slopes', max_slopes
#            print 'slope2', max_slopes2
#            print len(max_slopes)
#            print len(max_slopes2)
#            print delta
            #here calculate dzdt for each node, with initial time step
            #print "dwnst_nodes", dwnst_nodes
            for i in dwnst_nodes:
                if i==inlet_node:
                    qsin[i]=qsinlet
                    #print "qsin[inlet]", qsin[i]
                    #print "inlet area", drain_area[i]

                #calc deposition and erosion
                #dzver is vertical erosion/deposition only
#                print 'i ', i
#                print 'slope', max_slopes[i]
#                print 'area', drain_area[i]

                dep = alph*qsin[i]/drain_area[i]
                ero = -Kv * drain_area[i]**(0.5)*max_slopes[i]
#                print 'dep', dep
#                print 'ero', ero
                dzver[i] =  dep + ero
                #qsqt[i] = dep/-ero

                #calculate transport capacity
                qt[i]=Kv*drain_area[i]**(3./2.)*max_slopes[i]/alph
#                print 'drainarea', drain_area[i]
#                print 'maxslopes', max_slopes[i]
#                print 'ero', ero
#                print 'inode', i
                #lateral erosion component
                #potential lateral erosion initially set to 0
                petlat=0.                
                #bank height initially set to 0
                #z_bank=0.0
                #**********ADDED FOR WATER DEPTH CHANGE***************
                #may1, 2017. Need Q in m^3/s NOT M^3/yr!!!
                #water depth
#                qch=q[:]
#                qms=qch[i]/31536000.
#                wd=0.4*qms**0.35
                wd=0.4*(drain_area[i]*runoffms)**0.35
                if(0):
                    print 'i', i
                    print 'drain area[i]', drain_area[i]
                    print 'dr_area*runoffms', drain_area[i]*runoffms
                    print 'wd', wd
                
                #if node i flows downstream, continue. That is, if node i is the 
                #first cell at the top of the drainage network, don't go into this
                # loop because in this case, node i won't have a "donor" node found
                # in NodeFinder and needed to calculate the angle difference
                if Klr!= 0.0:
#                    print 'warning, in latero loop'
#                    print ' '
#                    print 'petlat before', petlat
#                    print 'i', i
                    if i in flowdirs:
                #if flowdirs[i] == i:
                    #Node_finder picks the lateral node to erode based on angle
                    # between segments between three nodes
                        [lat_node, inv_rad_curv]=Node_Finder2(grid, i, flowdirs, drain_area)
                    #node_finder returns the lateral node ID and the radius of curvature
                        lat_nodes[i]=lat_node
                        
                    #if the lateral node is not 0 continue. lateral node may be 
                    # 0 if a boundary node was chosen as a lateral node. then 
                    # radius of curavature is also 0 so there is no lateral erosion
                        if lat_node!=0.0:
                        #if the elevation of the lateral node is higher than primary node,
                        # calculate a new potential lateral erosion (L/T), which is negative
                            if z[lat_node] > z[i]:                           
                                petlat=-Kl*drain_area[i]*max_slopes[i]*inv_rad_curv
                            
                            #bank height. 
                            #z_bank=z[lat_node]-z[i]
                            
                            #the calculated potential lateral erosion is mutiplied by the length of the node 
                            #and the bank height, then added to an array, vol_lat_dt, for volume eroded 
                            #laterally  *per year* at each node. This vol_lat_dt is reset to zero for 
                            #each timestep loop. vol_lat_dt is added to itself more than one primary nodes are
                            # laterally eroding this lat_node                       
                                vol_lat_dt[lat_node]+=abs(petlat)*dx*wd                           
                        if (0):
                            print "i ", i
                            print 'lat_node', lat_node
                            print 'petlat after', petlat
                            print ' '
                # the following is always done, even if lat_node is 0 or lat_node 
                # lower than primary node. however, petlat is 0 in these cases
                    
                #send sediment downstream. sediment eroded from vertical incision 
                # and lateral erosion is sent downstream           	            	
                qsin[flowdirs[i]]+=qsin[i]-(dzver[i]*dx**2)-(petlat*dx*wd)   #qsin to next node
           
                #qsqt[i]=qsin[i]/qt[i]
                     

                if (0):
                    #if drain_area[i] > 500: i==inlet_node:
                    print 'node id', i
                    print "flowdirs[i]", flowdirs[i]
                    print "drain_area[i]", drain_area[i]
                    print "slope[i]", max_slopes[i]
                    print "qsin", qsin[i]
                    print "qsin downstream", qsin[flowdirs[i]]
                    print "dzver", dzver[i]
                    print 'petlat after', petlat
                    print 'wd', wd
                    print 'q[i]', q[i]
#                    print 'qms', qms
                    print " "
                    
            if(0):
            	print "dzlat", dzlat
                print "dzdt", dzdt
                
            dzdt=dzver            
            #Do a time-step check
            #If the downstream node is eroding at a slower rate than the 
            #upstream node, there is a possibility of flow direction reversal,
            #or at least a flattening of the landscape. 
            #Limit dt so that this flattening or reversal doesn't happen.
            #How close you allow these two points to get to eachother is
            #determined by the variable frac.
            
            #dtn is an arbitrarily large number to begin with, but will be adapted as we step through
            #the nodes
            dtn=dt*50 #starting minimum timestep for this round
            for i in dwnst_nodes:
                if(0):
                    print "node", i
                #are points converging? ie, downstream eroding slower than upstream
                dzdtdif = dzdt[flowdirs[i]]-dzdt[i]
                #if points converging, find time to zero slope
                if dzdtdif > 0. and max_slopes[i] > 1e-7:
                    dtflat = (z[i]-z[flowdirs[i]])/dzdtdif	#time to flat between points
                    #if time to flat is smaller than dt, take the lower value
                    if dtflat < dtn:
                        dtn = dtflat
                    #if dzdtdif*dtflat will make upstream lower than downstream, find time to flat
                    if dzdtdif*dtn > (z[i]-z[flowdirs[i]]):
                        dtn=(z[i]-z[flowdirs[i]])/dzdtdif
                        

            #print "out of ts loop"
            dtn*=frac			
            #new minimum timestep for this round of nodes
            dt=min(dtn, dt)
            #should now have a stable timestep.
#            print "stable time step=", dt
#            print delta
            
            #******Needed a stable timestep size to calculate volume eroded from lateral node for 
            # each stable time step********************
            
            #vol_lat is the total volume eroded from the lateral nodes through
            # the entire model run. So vol_lat is itself plus vol_lat_dt (for current loop)
            # times stable timestep size
            if(0):
                print "vol_lat before", vol_lat
                print "dt", dt
            vol_lat += vol_lat_dt*dt
            if (0): 
                print "vol_lat_dt", vol_lat_dt
                print "vol_lat after", vol_lat
                
            
            #this loop determines if enough lateral erosion has happened to change the height of the neighbor node.
            if Klr != 0.0:
#                print 'warning in lat ero, line 330'
                for i in dwnst_nodes:
                    lat_node=lat_nodes[i]
                    wd=0.4*(drain_area[i]*runoffms)**0.35
                    if lat_node!=0.0:
                        if z[lat_node] > z[i]:                        
                            
                            #September 11: changing so that voldiff is the volume that must be eroded 
                            # and the upper limit isn't the top of the node, but the water height at node i
                            # this would represent undercutting, slumping, and instant removal. 
                            #hmm, I have a mass balance problem here. 
                            voldiff=(z[i]+wd-z[flowdirs[i]])*dx**2
                            #if the total volume eroded from lat_node is greater than the volume 
                            # needed to be removed to make node equal elevation, 
                            # then instantaneously remove this height from lat node. already has timestep in it    
                            if vol_lat[lat_node]>=voldiff:
                                dzlat[lat_node]=z[flowdirs[i]]-z[lat_node]-0.001
                                if(0):
                                    print "chunk of lateral erosion occured"                        
                                    print "node", lat_node
                                    print "dzlat[lat_node]", dzlat[lat_node]
                                    print "vol_lat[latnode]", vol_lat[lat_node]
                                    print "vol_lat_dt", vol_lat_dt[lat_node]
                                    print "z[lat]", z[lat_node]
                                    print "z[i]", z[i]
                                    print "z[flowdirs[i]]", z[flowdirs[i]]
                                    print "volume diff", voldiff
                                    print 'wd', wd
                                    print "drain_area[i]", drain_area[i]
                                    print 'wd', wd
                                    print 'q[i]', q[i]
#                                    print 'qms', qms
#                                    print delta
                                #after the lateral node is eroded, reset its volume eroded to zero
                                vol_lat[lat_node]=0.0
            
            #multiply dzver(changed to dzdt above) by timestep size and combine with lateral erosion
            #dzlat, which is already a length for the chosen time step
            dz=dzdt*dt+dzlat        
            if(0):
                print 'dzlat[inlet_node]', dzlat[inlet_node]
                print 'dzdt[inlet_node]', dzdt[inlet_node]
                print 'dz[inlet_node]', dz[inlet_node]
            
            #change height of landscape
            z=dz+z
#            z[interior_nodes]=dz[interior_nodes]+z[interior_nodes]
            grid['node'][ 'topographic__elevation'] =  z
            #update elapsed time
            time=dt+time
#            print 'dz', dz
#            print 'time', time
                       
            #check to see that you are within 0.01% of the storm duration, if so done, if not continue

            if time > 0.9999*storm_dur:
                time = storm_dur
                #recalculate flow directions for plotting
#                flowdirs, drain_area, q, max_slopes, s, receiver_link = flow_router.route_flow(elevs=z, node_cell_area=node_area, runoff_rate=runoff)
                flow_router.route_flow(method='D8')
                s=grid.at_node['flow__upstream_node_order']
                drain_area_fr=grid.at_node['drainage_area']
                max_slopes=grid.at_node['topographic__steepest_slope']
                q=grid.at_node['surface_water__discharge']
                flowdirs=grid.at_node['flow__receiver_node']
                drain_area=q/dx**2    #this is the drainage area that I need for code below with an inlet set by spatially varible runoff.                
                #recalculate downstream order
                dsind = np.where((s >= min(interior_nodes)) & (s <= max(interior_nodes)))
                l=s[np.where((grid.status_at_node[s] == 0))[0]]
                dwnst_nodes=l
                dwnst_nodes=dwnst_nodes[::-1]
#                max_slopes=grid.calc_grad_at_active_link(z)
                #temporary hack for drainage area
#                drain_area[inlet_node]=inlet_area
            else:
                dt = storm_dur - time
                #recalculate flow directions
#                flowdirs, drain_area, q, max_slopes, s, receiver_link = flow_router.route_flow(elevs=z, node_cell_area=node_area, runoff_rate=runoff)
                flow_router.route_flow(method='D8')
                s=grid.at_node['flow__upstream_node_order']
                drain_area_fr=grid.at_node['drainage_area']
                max_slopes=grid.at_node['topographic__steepest_slope']
                q=grid.at_node['surface_water__discharge']
                flowdirs=grid.at_node['flow__receiver_node']                
                drain_area=q/dx**2    #this is the drainage area that I need for code below with an inlet set by spatially varible runoff.
                #recalculate downstream order
                dsind = np.where((s >= min(interior_nodes)) & (s <= max(interior_nodes)))
                l=s[np.where((grid.status_at_node[s] == 0))[0]]
                dwnst_nodes=l
                dwnst_nodes=dwnst_nodes[::-1]
                #what did I do here with this slopes, grad at active links?
#                max_slopes=grid.calc_grad_at_active_link(z)
                #temporary hack for drainage area
#                drain_area[inlet_node]=inlet_area
                #clear qsin for next loop
                qsin = grid.zeros(centering='node')
                qt = grid.zeros(centering='node')
                lat_nodes=np.zeros(grid.number_of_nodes)
                dzlat=np.zeros(grid.number_of_nodes)
                vol_lat_dt=np.zeros(grid.number_of_nodes)
                dzver=np.zeros(grid.number_of_nodes)
#        print 'dz=', dz[interior_nodes]
#        print 'dzdt= ', dzdt[interior_nodes]
#        print 'area', drain_area[interior_nodes]
#        print 'flowdirs', flowdirs[interior_nodes]
	        
    	
        return z, qt, qsin, dzdt, dzlat, flowdirs, drain_area, dwnst_nodes, max_slopes, dt       
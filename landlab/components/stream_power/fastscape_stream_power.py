#! /usr/env/python

"""
This module attempts to "component-ify" GT's Fastscape stream power erosion.
Created DEJH, March 2014.
"""

import numpy
from landlab import ModelParameterDictionary
UNDEFINED_INDEX = numpy.iinfo(int).max

class SPEroder(object):
    '''
    This class uses the Braun-Willett Fastscape approach to calculate the amount
    of erosion at each node in a grid, following a stream power framework.
    It needs to be supplied with the key variables:
        K_sp
        m_sp
        dt
    ...which it will draw from a supplied input file. n_sp has to be 1 for the
    BW algorithm to work.
        
    '''
    
    def __init__(self, grid, input_stream):
        self.grid = grid
        inputs = ModelParameterDictionary(input_stream)
        
        #User sets:
        self.K = inputs.read_float('K_sp')
        self.m = inputs.read_float('m_sp')
        try:
            self.n = inputs.read_float('n_sp')
        except:
            self.n = 1.
        try:
            self.dt = inputs.read_float('dt')
        except: #if dt isn't supplied, it must be set by another module, so look in the grid
            print 'Setting dynamic timestep from the grid...'
            self.dt = grid['timestep'] #this functionality doesn't exist yet, but it should
            
        #make storage variables
        self.A_to_the_m = grid.create_node_array_zeros()
        self.alpha = grid.empty(centering='node')
        
        if self.n != 1.:
            raise ValueError('The Braun Willett stream power algorithm requires n==1., sorry...')

    def erode(self, grid_in):
        self.grid = grid_in #the variables must be stored internally to the grid, in fields
        upstream_order_IDs = self.grid['node']['upstream_ID_order']
        #ordered_receivers = self.grid['node']['flow_receiver'][upstream_order_IDs]  #"j" in GT's sketch
        #nonboundaries = numpy.not_equal(upstream_order_IDs, ordered_receivers)
        z = self.grid['node']['planet_surface__elevation']
        #interior_nodes = numpy.greater_equal(self.grid['node']['links_to_flow_receiver'], -1)
        #interior_nodes = (self.grid['node']['links_to_flow_receiver'][upstream_order_IDs])[nonboundaries]
        #flow_link_lengths = self.grid.link_length[interior_nodes]
        flow_link_lengths = UNDEFINED_INDEX + numpy.zeros_like(self.alpha)
        defined_flow_receivers = numpy.greater_equal(self.grid['node']['links_to_flow_receiver'],-1)
        print self.grid.number_of_links
        print self.grid.number_of_nodes
        print (self.grid['node']['links_to_flow_receiver'])[defined_flow_receivers]
        flow_link_lengths[defined_flow_receivers] = self.grid.link_length[(self.grid['node']['links_to_flow_receiver'])[defined_flow_receivers]]
        numpy.power(self.grid['node']['drainage_area'], self.m, out=self.A_to_the_m)
        #self.alpha[nonboundaries] = self.K * self.dt * self.A_to_the_m[nonboundaries] / flow_link_lengths
        self.alpha = self.K * self.dt * self.A_to_the_m / flow_link_lengths


        for i in upstream_order_IDs:
            j = self.grid['node']['flow_receiver'][i]
            if i != j:
                z[i] = (z[i] + self.alpha[i]*z[j])/(1.0+self.alpha[i])
                        
#        for i in upstream_order_IDs[nonboundaries]: #this loop will be SLOW - is there a way of accelerating it using clever searching or cumsums?
#            j = self.grid['node']['flow_receiver'][i]   # receiver (downstream node) of i
#            print i,j
#            z[i] += self.alpha[i]*z[j]/(1.+self.alpha[i])
        
        self.grid['node']['planet_surface__elevation'] = z
        
        return self.grid


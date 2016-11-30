# -*- coding: utf-8 -*-
"""
Created on Mon Nov 28 08:52:08 2016

@author: Charlie
"""
import sys
import numpy as np
from landlab import Component
from landlab import FieldError

class HybridAlluvium(Component):
    """
    Stream Power with Alluvium Conservation and Entrainment (SPACE)
    
    Algorithm developed by G. Tucker, summer 2016.
    Component written by C. Shobe, 11/28/2016.
    """
    
    _name= 'HybridAlluvium'
    
    _input_var_names = (
        'flow__receiver_node',
        'flow__link_to_receiver_node',
        'flow__upstream_node_order',
        'topographic__steepest_slope',
        'drainage_area',
        'soil__depth'
    )
    
    _output_var_names = (
        'topographic__elevation'
        'soil__depth'
    )
    
    _var_units = {
        'flow__receiver_node': '-',
        'flow__link_to_receiver_node': '-',
        'flow__upstream_node_order': '-',
        'topographic__steepest_slope': '-',
        'drainage_area': 'm**2',
        'soil__depth': 'm',
        'topographic__elevation': 'm',
    }
    
    _var_mapping = {
        'flow__receiver_node': 'node',
        'flow__link_to_receiver_node': 'node',
        'flow__upstream_node_order': 'node',
        'topographic__steepest_slope': 'node',
        'drainage_area': 'node',
        'soil__depth': 'node',
        'topographic__elevation': 'node',
    }
    
    _var_doc = {
        'flow__receiver_node':
            'Node array of receivers (node that receives flow from current '
            'node)',
        'flow__link_to_receiver_node':
            'ID of link downstream of each node, which carries the discharge',
        'flow__upstream_node_order':
            'Node array containing downstream-to-upstream ordered list of '
            'node IDs',
        'topographic__steepest_slope':
            'Topographic slope at each node',
        'drainage_area':
            "Upstream accumulated surface area contributing to the node's "
            "discharge",
        'soil__depth': 
            'Depth of sediment above bedrock',
        'topographic__elevation': 
            'Land surface topographic elevation',
    }
    
    def __init__(self, grid, K_sed=None, K_br=None, a=None, 
                 F_f=None, phi=None, H_star=None, v_s=None, **kwds):
        """Initialize the HybridAlluvium model.
        
        Parameters
        ----------
        grid : ModelGrid
            Landlab ModelGrid object
        K_sed : float
            Erodibility constant for sediment (units vary).
        K_br : float
            Erodibility constant for bedrock (units vary).
        a : float
            Exponent on drainage area to give discharge.
        F_f : float
            Fraction of permanently suspendable fines in bedrock [-].
        phi : float
            Sediment porosity [-].
        H_star : float
            Sediment thickness required for full entrainment [L].
        v_s : float
            Settling velocity for chosen grain size metric [L/T].
        """
        
        if 'soil__depth' not in grid.at_node:        
            grid.add_zeros('soil__depth', at='node', dtype=float)
        self._grid = grid #store grid
        self.K_sed = float(K_sed)
        self.K_br = float(K_br)
        self.a = float(a)
        self.F_f = float(F_f)
        self.phi = float(phi)
        self.H_star = float(H_star)
        self.v_s = float(v_s)
        self.drainage_area = grid.at_node['drainage_area']
        self.flow_receivers = grid.at_node['flow__receiver_node']
        self.stack = grid.at_node['flow__upstream_node_order']
        self.elevation = grid.at_node['topographic__elevation']
        self.slope = grid.at_node['topographic__steepest_slope']
        self.soil__depth = grid.at_node['soil__depth']
        self.stack_links = grid.at_node['flow__link_to_receiver_node']
        self.link_lengths = grid._length_of_link_with_diagonals
        self.upstream_node_order = grid.at_node['flow__upstream_node_order']
        self.bedrock__elevation = grid.add_zeros(
            'bedrock__elevation', at='node', dtype=float)
        self.qs = np.zeros(self.grid.number_of_nodes)#take out of dt
        
    def run_one_step(self, dt=1.0, **kwds):
        """Calculate change in rock and alluvium thickness for
           a time period 'dt'.
        """
        self.q = (self.drainage_area ** self.a) / self.grid.width_of_face[0] #MEANS ONLY WORKS ON RASTER GRID FOR NOW
        
        self.Es = self.K_sed * self.q * self.slope * \
            (1.0 - np.exp(-self.soil__depth / self.H_star)) #erosion of sediment at every node, vectorized
        if sum(np.isnan(self.Es)):
            sys.exit('first NaN occurrence 1')
        self.Er = self.K_br * self.q * self.slope * np.exp(-self.soil__depth / self.H_star)
        if sum(np.isnan(self.Er)):
            sys.exit('first NaN occurrence 2')
        #self.qs = np.zeros(self.grid.number_of_nodes)#take out of dt
        if sum(np.isnan(self.qs)):
            sys.exit('first NaN occurrence 3')
        self.qs_in = np.zeros(self.grid.number_of_nodes)
        if sum(np.isnan(self.qs_in)):
            sys.exit('first NaN occurrence 4')
        for j in np.flipud(self.stack):
            if self.q[j] == 0:
                self.qs[j] = 0
            else:
                self.qs[j] = (((self.Es[j]) + (1-self.F_f) * self.Er[j]) / (self.v_s / self.q[j])) * \
                    (1.0 - np.exp(-self.link_lengths[j] * self.v_s / self.q[j])) + (self.qs_in[j] * \
                    np.exp(-self.link_lengths[j] * self.v_s / self.q[j]))
            self.qs_in[self.flow_receivers[j]] += self.qs[j]
        self.soil__depth[self.soil__depth < 1e-9] = 0
        self.soil__depth[self.q > 0] += dt * ((self.qs[self.q > 0] * (self.v_s / self.q[self.q > 0])) - self.K_sed * self.q[self.q > 0] * \
            self.slope[self.q > 0] *(1.0 - np.exp(-self.soil__depth[self.q > 0]/self.H_star)))
        print 'depth 8198', self.soil__depth[8198]#min(self.soil__depth)
        print max(self.soil__depth)
        if np.any(self.soil__depth < 0):
            sys.exit('negative soil')
        self.bedrock__elevation[self.q > 0] += dt * (self.K_br * self.q[self.q > 0] * self.slope[self.q > 0] * (np.exp(-self.soil__depth[self.q > 0] /\
            self.H_star))) #took out uplift b/c that is user's job outside this component
            
        self.elevation = self.bedrock__elevation + self.soil__depth 
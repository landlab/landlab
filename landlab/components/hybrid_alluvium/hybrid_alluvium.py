# -*- coding: utf-8 -*-
"""
Created on Mon Nov 28 08:52:08 2016

@author: Charlie
"""
import numpy as np
from landlab import Component

class HybridAlluvium(Component):
    """
    Stream Power with Alluvium Conservation and Entrainment (SPACE)
    
    Algorithm developed by G. Tucker, summer 2016.
    Component written by C. Shobe, begun 11/28/2016.
    
    Currently only works on square raster grids (3/14/2017).
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
                 F_f=None, phi=None, H_star=None, v_s=None, 
                 m_sp=None, n_sp=None, **kwds):
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
        if 'water__discharge' not in grid.at_node:        
            grid.add_zeros('water__discharge', at='node', dtype=float)
        self._grid = grid #store grid
        self.m_sp = float(m_sp)
        self.n_sp = float(n_sp)
        self.K_sed = float(K_sed)
        self.K_br = float(K_br)
        self.a = float(a)
        self.F_f = float(F_f)
        self.phi = float(phi)
        self.H_star = float(H_star)
        self.v_s = float(v_s)
        self.drainage_area = grid.at_node['drainage_area']
        self.q = grid.at_node['water__discharge']
        self.flow_receivers = grid.at_node['flow__receiver_node']
        self.stack = grid.at_node['flow__upstream_node_order']
        self.topographic__elevation = grid.at_node['topographic__elevation']
        self.slope = grid.at_node['topographic__steepest_slope']
        self.soil__depth = grid.at_node['soil__depth']
        self.stack_links = grid.at_node['flow__link_to_receiver_node']
        self.link_lengths = grid._length_of_link_with_diagonals
        self.upstream_node_order = grid.at_node['flow__upstream_node_order']
        self.bedrock__elevation = grid.add_zeros(
            'bedrock__elevation', at='node', dtype=float)
        self.bedrock__elevation += self.topographic__elevation.copy()
        self.qs = grid.add_zeros('sediment__flux', at='node', dtype=float)

    #def simple_shear_stress(self):
    def simple_stream_power(self):
        self.Es = self.K_sed * np.power(self.grid.at_node['drainage_area'], 
                                        self.m_sp) * \
            np.power(self.slope, self.n_sp) * \
            (1.0 - np.exp(-self.soil__depth / self.H_star))
        self.Er = self.K_br * np.power(self.grid.at_node['drainage_area'], 
                                       self.m_sp) *\
            np.power(self.slope, self.n_sp) * \
            np.exp(-self.soil__depth / self.H_star)
        self.q = np.power(self.grid.at_node['drainage_area'], self.m_sp)
        self.sed_erosion_term = self.K_sed * \
            np.power(self.grid.at_node['drainage_area'], self.m_sp) * \
            np.power(self.slope, self.n_sp)
        self.br_erosion_term = self.K_br * \
            np.power(self.grid.at_node['drainage_area'], self.m_sp) * \
            np.power(self.slope, self.n_sp)
        return self.Es
    #def threshold_shear_stress(self, ):
        
    def run_one_step(self, dt=1.0, flooded_nodes=None, method=None, **kwds):
        """Calculate change in rock and alluvium thickness for
           a time period 'dt'.
        """
        dt = dt * 3.154e7 #so input dt is in years
        
        #Choose a method for calculating erosion:
        if method == 'simple_shear_stress':        
            self.simple_shear_stress()        
        elif method == 'simple_stream_power':
            self.simple_stream_power()
        elif method == 'threshold_shear_stress':
            self.threshold_shear_stress()
        else:
            raise ValueError('Specify an erosion method!')
            
        self.qs_in = np.zeros(self.grid.number_of_nodes)            
            
        for j in np.flipud(self.stack):
            if self.q[j] == 0:
                self.qs[j] = 0
            else:
                self.qs[j] = (((self.Es[j]) + (1-self.F_f) * self.Er[j]) / \
                    (self.v_s / self.q[j])) * (1.0 - \
                    np.exp(-self.link_lengths[j] * self.v_s / self.q[j])) + \
                    (self.qs_in[j] * np.exp(-self.link_lengths[j] * \
                    self.v_s / self.q[j]))
            self.qs_in[self.flow_receivers[j]] += self.qs[j]
        erosion_pertime = np.zeros(self.grid.number_of_nodes)
        erosion = np.zeros(self.grid.number_of_nodes)
        deposition_pertime = np.zeros(self.grid.number_of_nodes)
        deposition = np.zeros(self.grid.number_of_nodes)
        erosion_pertime = self.simple_stream_power()
        erosion[self.q > 0] = dt * erosion_pertime[self.q > 0]
        deposition_pertime[self.q > 0] = (self.qs[self.q > 0] * \
            (self.v_s / self.q[self.q > 0]))
        deposition[self.q > 0] = dt * deposition_pertime[self.q > 0]
        
        #now, the analytical solution to soil thickness in time:
        #need to distinguish D=kqS from all other cases to save from blowup!
        
        flooded = self._grid.nodes.flatten() == flooded_nodes
        
        #distinguish cases:
        blowup = deposition_pertime == self.K_sed * self.q * self.slope

        ##first, potential blowup case:
        #positive slopes, not flooded
        self.soil__depth[(self.q > 0) & (blowup==True) & (self.slope > 0) & \
            (flooded==False)] = self.H_star * \
            np.log((erosion_pertime[(self.q > 0) & (blowup==True) & \
            (self.slope > 0) & (flooded==False)] / self.H_star) * dt + \
            np.exp(self.soil__depth[(self.q > 0) & (blowup==True) & \
            (self.slope > 0) & (flooded==False)] / self.H_star))
        #positive slopes, flooded
        self.soil__depth[(self.q > 0) & (blowup==True) & (self.slope > 0) & \
            (flooded==True)] = (deposition_pertime[(self.q > 0) & \
            (blowup==True) & (flooded==True)] / (1 - self.phi)) * dt   
        #non-positive slopes, not flooded
        self.soil__depth[(self.q > 0) & (blowup==True) & (self.slope <= 0) & \
            (flooded==False)] += (deposition_pertime[(self.q > 0) & \
            (blowup==True) & (self.slope <= 0) & (flooded==False)] / \
            (1 - self.phi)) * dt    
        
        ##more general case:
        #positive slopes, not flooded
        self.soil__depth[(self.q > 0) & (blowup==False) & (self.slope > 0) & \
            (flooded==False)] = self.H_star * \
            np.log((1 / ((deposition_pertime[(self.q > 0) & (blowup==False) & \
            (self.slope > 0) & (flooded==False)] * (1 - self.phi)) / \
            (self.sed_erosion_term[(self.q > 0) & (blowup==False) & \
            (self.slope > 0) & (flooded==False)]) - 1)) * \
            (np.exp((deposition_pertime[(self.q > 0) & (blowup==False) & \
            (self.slope > 0) & (flooded==False)] * (1 - self.phi) - \
            (self.sed_erosion_term[(self.q > 0) & (blowup==False) & \
            (self.slope > 0) & (flooded==False)]))*(dt / self.H_star)) * \
            (((deposition_pertime[(self.q > 0) & (blowup==False) & \
            (self.slope > 0) & (flooded==False)] * (1 - self.phi) / \
            (self.sed_erosion_term[(self.q > 0) & (blowup==False) & \
            (self.slope > 0) & (flooded==False)])) - 1) * \
            np.exp(self.soil__depth[(self.q > 0) & (blowup==False) & \
            (self.slope > 0) & (flooded==False)] / self.H_star)  + 1) - 1))
        #places where slope <= 0 but not flooded:
        self.soil__depth[(self.q > 0) & (blowup==False) & (self.slope <= 0) & \
            (flooded==False)] += (deposition_pertime[(self.q > 0) & \
            (blowup==False) & (self.slope <= 0) & (flooded==False)] / \
            (1 - self.phi)) * dt     
        #flooded nodes:        
        self.soil__depth[(self.q > 0) & (blowup==False) & (flooded==True)] += \
            (deposition_pertime[(self.q > 0) & (blowup==False) & \
            (flooded==True)] / (1 - self.phi)) * dt     

        #check against negative soil thickness
        self.soil__depth[self.soil__depth + deposition - erosion < 0] = 0
        
        if np.any(self.soil__depth < 0):
            raise ValueError('Negative soil! Killing model.')
        self.bedrock__elevation[self.q > 0] += dt * \
            (-self.br_erosion_term[self.q > 0] * \
            (np.exp(-self.soil__depth[self.q > 0] / self.H_star)))

        self.topographic__elevation = self.bedrock__elevation + \
            self.soil__depth 
        
        #save as grid fields
        self._grid['node']['topographic__elevation'] = \
            self.topographic__elevation
        self._grid['node']['bedrock__elevation'] = self.bedrock__elevation
        self._grid['node']['soil__depth'] = self.soil__depth
        self._grid['node']['water__discharge'] = self.q
        self._grid['node']['sediment__flux'] = self.qs
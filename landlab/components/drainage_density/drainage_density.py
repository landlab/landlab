# -*- coding: utf-8 -*-
"""
Created on Mon Jul 11 10:00:25 2016

@author: Charlie Shobe

for West Valley Project: Landlab component to calculate drainage density

"""
from landlab import Component
import numpy as np

class DrainageDensity(Component):
    """
    Calculate drainage density overa DEM.
    
    Landlab component that implements the distance to channel algorithm of
    Tucker et al., 2001.
    
    Written by C. Shobe on 7/11/2016
    
    Construction:
    
        DrainageDensity(grid, channel_network_name='string')
        
    Parameters
    ----------
    grid : ModelGrid
    channel_network_name : string naming a boolean grid field with 'True'
        in channels and 'False' elsewhere.
        
    Examples
    ---------
    >>> import numpy as np
    >>> from landlab import RasterModelGrid
    >>> from landlab.components.flow_routing import FlowRouter
    >>> from landlab.components import FastscapeEroder
    >>> mg = RasterModelGrid((10, 10), 1.0)
    >>> _ = mg.add_zeros('node', 'topographic__elevation')
    >>> np.random.seed(50)
    >>> noise = np.random.rand(100)
    >>> mg['node']['topographic__elevation'] += noise
    >>> mg['node']['topographic__elevation'] # doctest: +NORMALIZE_WHITESPACE
    array([ 0.49460165,  0.2280831 ,  0.25547392,  0.39632991,  0.3773151 ,
        0.99657423,  0.4081972 ,  0.77189399,  0.76053669,  0.31000935,
        0.3465412 ,  0.35176482,  0.14546686,  0.97266468,  0.90917844,
        0.5599571 ,  0.31359075,  0.88820004,  0.67457307,  0.39108745,
        0.50718412,  0.5241035 ,  0.92800093,  0.57137307,  0.66833757,
        0.05225869,  0.3270573 ,  0.05640164,  0.17982769,  0.92593317,
        0.93801522,  0.71409271,  0.73268761,  0.46174768,  0.93132927,
        0.40642024,  0.68320577,  0.64991587,  0.59876518,  0.22203939,
        0.68235717,  0.8780563 ,  0.79671726,  0.43200225,  0.91787822,
        0.78183368,  0.72575028,  0.12485469,  0.91630845,  0.38771099,
        0.29492955,  0.61673141,  0.46784623,  0.25533891,  0.83899589,
        0.1786192 ,  0.22711417,  0.65987645,  0.47911625,  0.07344734,
        0.13896007,  0.11230718,  0.47778497,  0.54029623,  0.95807105,
        0.58379231,  0.52666409,  0.92226269,  0.91925702,  0.25200886,
        0.68263261,  0.96427612,  0.22696165,  0.7160172 ,  0.79776011,
        0.9367512 ,  0.8537225 ,  0.42154581,  0.00543987,  0.03486533,
        0.01390537,  0.58890993,  0.3829931 ,  0.11481895,  0.86445401,
        0.82165703,  0.73749168,  0.84034417,  0.4015291 ,  0.74862   ,
        0.55962945,  0.61323757,  0.29810165,  0.60237917,  0.42567684,
        0.53854438,  0.48672986,  0.49989164,  0.91745948,  0.26287702])
    >>> fr = FlowRouter(mg)
    >>> fsc = FastscapeEroder(mg, K_sp=.01, m_sp=.5, n_sp=1)
    >>> for x in range(100):
    ...     fr.run_one_step()
    ...     fsc.run_one_step(dt = 10.0)
    ...     mg.at_node['topographic__elevation'][mg.core_nodes] += .01
    >>> channels = mg['node']['drainage_area'] > 5
    >>> _ = mg.add_field('node', 'channel_network', channels)
    >>> dd = DrainageDensity(mg, channel_network_name='channel_network')
    >>> mean_drainage_density = dd.calc_drainage_density()
    >>> print np.around(mean_drainage_density, 10)
    0.3831100571
    """
    
    _name = 'DrainageDensity'
    
    _input_var_names = (
        'topographic__elevation'
        'topographic__gradient'
        'drainage_area'
    )
    
    _output_var_names = (
        'distance_to_channel'
        'drainage_density' #this is the CSDMS standard name
    )
    
    _var_units = {
        'topographic__elevation': 'm',
        'topographic__gradient': 'm/m',
        'drainage_area': 'm^2',
        'distance_to_channel': 'm',
        'drainage_density': '1/m'
    }
    
    _var_mapping = {
        'topographic__elevation': 'node',
        'topographic__gradient': 'link',
        'drainage_area': 'node',
        'distance_to_channel': 'node',
        'drainage_density': 'node'
    }
    
    _var_doc = {
        'topographic__elevation':
            'elevation of the ground surface relative to some datum',
        'topographic__gradient':
            'gradient of the ground surface',
        'drainage_area':
            'contributing drainage area',
        'distance_to_channel':
            'distance from each node to the nearest channel',
        'drainage_density':
            'total length of channels per unit area'
    }
    
    def __init__(self, grid, channel_network_name = 'channel_network', **kwds):
        """Initialize the DrainageDensity component.
        
        Parameters
        ----------
        grid : ModelGrid
            Landlab ModelGrid object
        channel_network_name : string, optional (defaults to 'channel_network')
            Name of Landlab grid field that holds 1's where channels exist and 0's elsewhere
        """
        
        #Store grid
        self._grid = grid
        
        # Create fields...
        
        #   Channel network
        if channel_network_name in grid.at_node:
            self.channel_network = grid.at_node['channel_network']
        else:
            self.channel_network = grid.add_zeros('node', 'channel_network')

        #   Flow receivers
        if 'flow__receiver_node' in grid.at_node:
            self.flow_receivers = grid.at_node['flow__receiver_node']
        else:
            self.flow_receivers = grid.add_zeros('node', 'flow__receiver_node')
            
        #   Links to receiver nodes
        if 'flow__link_to_receiver_node' in grid.at_node: #flow__link_to_receiver_node
            self.stack_links = grid.at_node['flow__link_to_receiver_node']
        else:
            self.stack_links = grid.add_zeros('node', 'flow__link_to_receiver_node')
        
        #   Distance to channel
        if 'distance_to_channel' in grid.at_node:
            self.distance_to_channel = grid.at_node['distance_to_channel']
        else:
            self.distance_to_channel = grid.add_zeros('node', 'distance_to_channel')
        
    def calc_drainage_density(self, **kwds): #there is no 'run_one_step' methid b/c this is a tool, not a model.
        """Calculate distance to channel and drainage density, after
        Tucker et al., 2001.
        """
        for node in range(self.grid.number_of_nodes):
            if self.channel_network[node] == True: #we're standing in a channel
                distance = 0
            else:
                distance = 0
                flag = 0
                node_iter = node
                while flag == 0: #as long as we haven't hit a channel yet...
                    if self.flow_receivers[node_iter] == node_iter: #if no flow receiver (boundary probably)
                        self.channel_network[node_iter] = True #convince the node it's a channel
                    else:
                        distance += self.grid._length_of_link_with_diagonals[self.stack_links[node_iter]]
                        node_iter = self.flow_receivers[node_iter]
                    if self.channel_network[node_iter] == True: #we've hit a channel
                        self.distance_to_channel[node] = distance #save distance to channel
                        flag = 1
        self.grid['node']['distance_to_channel'] = self.distance_to_channel
        landscape_drainage_density = 1 / (2.0 * np.mean(self.grid.at_node['distance_to_channel'][self.grid.core_nodes]))#self.distance_to_channel)) #this is THE drainage density
        return landscape_drainage_density
                
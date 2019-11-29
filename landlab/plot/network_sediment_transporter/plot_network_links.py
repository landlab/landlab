# -*- coding: utf-8 -*-
"""
Created on Fri Oct 18 11:57:55 2019

This code plots the network and colors each link according to a link attribute.

@author: jczuba
"""
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.collections import LineCollection
    
def plot_network_links(grid, link_attribute, *args, **kwargs):
    # Need to verify that the link_attribute is valid, 
    # error messages if that's not the case
    
    if np.shape(link_attribute) != np.shape(grid.at_link["channel_slope"]):
        msg = "Plot Network Links: attribute must be of shape (links,)"
        raise ValueError(msg)
    
    # %% plot attribute as colored link of network
    
    # We need to set the plot limits, they will not autoscale
    fig, ax = plt.subplots()
    ax.set_xlim(grid.x_of_node.min(), grid.x_of_node.max())
    ax.set_ylim(grid.y_of_node.min(), grid.y_of_node.max())
       
#    line_segments = LineCollection(segments)
#    line_segments.set_array(link_attribute)
#    ax.add_collection(line_segments)
#    axcb = fig.colorbar(line_segments)
#        
#    plt.show()
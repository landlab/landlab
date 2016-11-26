# -*- coding: utf-8 -*-
"""
stream_power_smooth_threshold.py: Defines the StreamPowerSmoothThresholdEroder,
which is a version of the FastscapeEroder (derived from it).

StreamPowerSmoothThresholdEroder uses a mathematically smooth threshold 
formulation, rather than one with a singularity. The erosion rate is defined as
follows:

$\omega = K A^m S$

$E = \omega - \omega_c \left[ 1 - \exp ( -\omega / \omega_c ) \right]$

Created on Sat Nov 26 08:36:49 2016

@author: gtucker
"""

from landlab.components import FastscapeEroder
import numpy as np

UNDEFINED_INDEX = -1



class StreamPowerSmoothThresholdEroder(FastscapeEroder):
    """Stream erosion component with smooth threshold function.
    
    [ADD DOCS AND TESTS HERE]
    """
    
    def __init__(self, grid, K_sp=None, m_sp=0.5, n_sp=1., threshold_sp=1.,
                 rainfall_intensity=1., **kwargs):
        """Initialize StreamPowerSmoothThresholdEroder."""
        
        # Call base-class init
        super(StreamPowerSmoothThresholdEroder, 
              self).__init__(grid, K_sp=None, m_sp=0.5, n_sp=1.,
                             threshold_sp=1., rainfall_intensity=1., **kwargs)

    def run_one_step(self, dt, flooded_nodes=None,
                     rainfall_intensity_if_used=None, **kwds):
        """Run one forward iteration of duration dt."""
        
        # Set up needed arrays
        #
        # Get shorthand for elevation field ("z"), and for up-to-downstream
        # ordered array of node IDs ("upstream_order_IDs")
        upstream_order_IDs = self._grid['node']['flow__upstream_node_order']
        z = self._grid['node']['topographic__elevation']

        # Get an array of flow-link length for each node that has a defined
        # receiver (i.e., that drains to another node).
        defined_flow_receivers = np.not_equal(self._grid['node'][
            'flow__link_to_receiver_node'], UNDEFINED_INDEX)
        flow_link_lengths = self._grid._length_of_link_with_diagonals[
            self._grid['node']['flow__link_to_receiver_node'][
                defined_flow_receivers]]
                
        # (Later on, add functionality for a runoff rate, or discharge, or
        # variable K)
                
        # Set up alpha, beta, delta arrays

        # Iterate over nodes from downstream to upstream, using scipy's
        # 'newton' function to find new elevation at each node in turn.
        
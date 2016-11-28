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

if __name__ == '__main__':
    from landlab.components import FastscapeEroder
else:
    from .fastscape_stream_power import FastscapeEroder
import numpy as np
from scipy.optimize import newton

UNDEFINED_INDEX = -1




def new_elev(x, a, b, c, d):
    """Equation for elevation of a node at timestep t+1.
    
    Parameters
    ----------
    x : float
        Value of new elevation
    a : float
        Parameter = K A^m dt / L
    b : float
        Parameter = K A^m / wc L
    c : float
        Elevation of downstream node, z_j
    d : float
        z(t) + a z_j + wc
    """
    return x * (1.0 + a) + np.exp(-b * (x - c)) - d


def new_elev_prime(x, a, b, c, d=0.0):
    """Equation for elevation of a node at timestep t+1.
    
    Parameters
    ----------
    x : float
        Value of new elevation
    a : float
        Parameter = K A^m dt / L
    b : float
        Parameter = K A^m / wc L
    c : float
        Elevation of downstream node, z_j
    d : n/a
        Placeholder; not used
    """
    return (1.0 + a) - b * np.exp(-b * (x - c))


def new_elev_prime2(x, a, b, c, d=0.0):
    """Equation for elevation of a node at timestep t+1.
    
    Parameters
    ----------
    x : float
        Value of new elevation
    a : n/a
        Placeholder; not used
    b : float
        Parameter = K A^m / wc L
    c : float
        Elevation of downstream node, z_j
    d : n/a
        Placeholder; not used
    """
    return b * b * np.exp(-b * (x - c))



class StreamPowerSmoothThresholdEroder(FastscapeEroder):
    """Stream erosion component with smooth threshold function.
    
    Parameters
    ----------
    
    
    Examples
    --------
    >>> from landlab import RasterModelGrid
    >>> rg = RasterModelGrid((3, 3), 1.0)
    >>> z = rg.add_zeros('node', 'topographic__elevation')
    >>> from landlab.components import StreamPowerSmoothThresholdEroder
    >>> sp = StreamPowerSmoothThresholdEroder(rg, K_sp=1.0)
    >>> sp.thresholds
    1.0
    """
    
    def __init__(self, grid, K_sp=None, m_sp=0.5, n_sp=1., threshold_sp=1.,
                 rainfall_intensity=1., **kwargs):
        """Initialize StreamPowerSmoothThresholdEroder."""
        
        # Call base-class init
        super(StreamPowerSmoothThresholdEroder, 
              self).__init__(grid, K_sp=K_sp, m_sp=m_sp, n_sp=n_sp, 
                             threshold_sp=threshold_sp,
                             rainfall_intensity=rainfall_intensity, **kwargs)

        # Arrays with parameters for use in implicit solver
        self.beta = grid.empty(at='node')
        self.delta = grid.empty(at='node')


    def run_one_step(self, dt, flooded_nodes=None, runoff_rate=None, **kwds):
        """Run one forward iteration of duration dt.
        
        Parameters
        ----------
        dt : float
            Time step size
        flooded_nodes : ndarray of int (optional)
            Indices of nodes in lakes/depressions
        runoff_rate : (not used yet)
            (to be added later)
            
        Examples
        --------
        >>> from landlab import RasterModelGrid
        >>> rg = RasterModelGrid((3, 3), 1.0)
        >>> rg.set_closed_boundaries_at_grid_edges(False, True, True, True)
        >>> z = rg.add_zeros('node', 'topographic__elevation')
        >>> z[4] = 1.0
        >>> from landlab.components import FlowRouter
        >>> fr = FlowRouter(rg, method='D4')
        >>> fr.run_one_step()
        >>> from landlab.components import StreamPowerSmoothThresholdEroder
        >>> sp = StreamPowerSmoothThresholdEroder(rg, K_sp=1.0)
        >>> sp.run_one_step(dt=1.0)
        >>> sp.alpha
        array([ 0.,  0.,  0.,  0.,  1.,  0.,  0.,  0.,  0.])
        >>> sp.beta
        array([ 0.,  0.,  0.,  0.,  1.,  0.,  0.,  0.,  0.])
        >>> sp.gamma
        array([ 0.,  0.,  0.,  0.,  2.,  0.,  0.,  0.,  0.])
        """
        
        # Set up needed arrays
        #
        # Get shorthand for elevation field ("z"), and for up-to-downstream
        # ordered array of node IDs ("upstream_order_IDs")
        upstream_order_IDs = self._grid['node']['flow__upstream_node_order']
        z = self._grid['node']['topographic__elevation']
        flow_receivers = self._grid['node']['flow__receiver_node']

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
        #
        #   First, compute drainage area raised to the power m.
        np.power(self._grid['node']['drainage_area'], self.m,
                 out=self.A_to_the_m)
                 
        #   Alpha
        self.alpha[defined_flow_receivers] = (self.K * dt
            * self.A_to_the_m[defined_flow_receivers] / flow_link_lengths)

        #   Beta
        self.beta[defined_flow_receivers] = ((self.K
            * self.A_to_the_m[defined_flow_receivers] ) 
            / (self.thresholds * flow_link_lengths))

        #   Gamma
        self.delta[defined_flow_receivers] = (z[defined_flow_receivers]
            + (self.alpha[defined_flow_receivers]
               * z[flow_receivers[defined_flow_receivers]]) 
            + self.thresholds)

        # Iterate over nodes from downstream to upstream, using scipy's
        # 'newton' function to find new elevation at each node in turn.
        for node in upstream_order_IDs:
            print(node)
            if defined_flow_receivers[node]:
                z[node] = newton(new_elev, z[node],
                                 fprime=new_elev_prime,
                                 args=(self.alpha[node], self.beta[node],
                                       z[flow_receivers[node]],
                                       self.delta[node]), 
                                 fprime2=new_elev_prime2)
                print(node, z[node])

        # TODO: handle case self.thresholds = 0




if __name__ == '__main__':
    from landlab import RasterModelGrid
    rg = RasterModelGrid((3, 4), 1.0)
    z = rg.add_zeros('node', 'topographic__elevation')
    rg.set_closed_boundaries_at_grid_edges(False, True, True, True)
    z[5] = 2.0
    z[6] = 1.0
    from landlab.components import FlowRouter
    fr = FlowRouter(rg, method='D4')
    fr.run_one_step()
    sp = StreamPowerSmoothThresholdEroder(rg, K_sp=1.0)
    sp.run_one_step(dt=1.0)

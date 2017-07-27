#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Thu Jul 27 14:23:25 2017

@author: gtucker
"""

from landlab import RasterModelGrid
from landlab.components import ErosionDeposition, FlowAccumulator
import numpy as np
from numpy.testing import assert_equal

def test_erodep_slope_area():
    """Test steady state run with Vs << 1, Vs >> 1, and Vs = 1."""

    # Set up a 5x5 grid with open boundaries and low initial elevations.
    rg = RasterModelGrid((5, 5))
    z = rg.add_zeros('node', 'topographic__elevation')
    z[:] = 0.01 * rg.x_of_node

    # Create a D8 flow handler
    fa = FlowAccumulator(rg, flow_director='FlowDirectorD8')

    # Parameter values for test 1
    K = 0.001
    vs = 0.0001
    U = 0.001
    dt = 10.0

    # Create the ErosionDeposition component...
    ed = ErosionDeposition(rg, K=K, phi=0.0, v_s=vs, m_sp=0.5, n_sp=1.0,
                           method='simple_stream_power',
                           discharge_method='drainage_area', 
                           area_field='drainage_area',
                           solver='adaptive')

    # ... and run it to steady state.
    for i in range(1000):
        fa.run_one_step()
        ed.run_one_step(dt=dt)
        z[rg.core_nodes] += U * dt

    # Test the results
    s = rg.at_node['topographic__steepest_slope']
    sa_factor = (1.0 + vs) * U / K
    a11 = 2.0
    a12 = 1.0
    s = rg.at_node['topographic__steepest_slope']    
    s11 = sa_factor * (a11 ** -0.5)
    s12 = sa_factor * (a12 ** -0.5)
    assert_equal(np.round(s[11], 3), np.round(s11, 3))
    assert_equal(np.round(s[12], 3), np.round(s12, 3))

    # Next test: big Vs
    K = 1.0
    vs = 1000.0
    
    # Create the ErosionDeposition component...
    ed = ErosionDeposition(rg, K=K, phi=0.0, v_s=vs, m_sp=0.5, n_sp=1.0,
                           method='simple_stream_power',
                           discharge_method='drainage_area', 
                           area_field='drainage_area',
                           solver='adaptive')

    # ... and run it to steady state.
    for i in range(1000):
        fa.run_one_step()
        ed.run_one_step(dt=dt)
        z[rg.core_nodes] += U * dt

    # Test the results
    sa_factor = (1.0 + vs) * U / K
    s11 = sa_factor * (a11 ** -0.5)
    s12 = sa_factor * (a12 ** -0.5)
    assert_equal(np.round(s[11], 3), np.round(s11, 3))
    assert_equal(np.round(s[12], 3), np.round(s12, 3))

    # Final test: Vs = 1
    K = 0.002
    vs = 1.0

    # Create the ErosionDeposition component...
    ed = ErosionDeposition(rg, K=K, phi=0.0, v_s=vs, m_sp=0.5, n_sp=1.0,
                           method='simple_stream_power',
                           discharge_method='drainage_area', 
                           area_field='drainage_area',
                           solver='adaptive')

    # ... and run it to steady state.
    for i in range(1000):
        fa.run_one_step()
        ed.run_one_step(dt=dt)
        z[rg.core_nodes] += U * dt

    # Test the results
    sa_factor = (1.0 + vs) * U / K
    s11 = sa_factor * (a11 ** -0.5)
    s12 = sa_factor * (a12 ** -0.5)
    assert_equal(np.round(s[11], 3), np.round(s11, 3))
    assert_equal(np.round(s[12], 3), np.round(s12, 3))


if __name__ == '__main__':
    test_erodep_slope_area()

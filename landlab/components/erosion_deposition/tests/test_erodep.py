#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Thu Jul 27 14:23:25 2017

@author: gtucker
"""

from landlab import RasterModelGrid, HexModelGrid
from landlab.components import ErosionDeposition, FlowAccumulator
import numpy as np
from numpy.testing import assert_equal

def test_can_run_with_hex():
    """Test that model can run with hex model grid."""

    # Set up a 5x5 grid with open boundaries and low initial elevations.
    mg = HexModelGrid(7, 7)
    z = mg.add_zeros('node', 'topographic__elevation')
    z[:] = 0.01 * mg.x_of_node

    # Create a D8 flow handler
    fa = FlowAccumulator(mg, flow_director='FlowDirectorSteepest')

    # Parameter values for test 1
    K = 0.001
    vs = 0.0001
    U = 0.001
    dt = 10.0

    # Create the ErosionDeposition component...
    ed = ErosionDeposition(mg, K=K, phi=0.0, v_s=vs, m_sp=0.5, n_sp=1.0,
                           method='simple_stream_power',
                           discharge_method='drainage_area',
                           area_field='drainage_area',
                           solver='adaptive')

    # ... and run it to steady state.
    for i in range(2000):
        fa.run_one_step()
        ed.run_one_step(dt=dt)
        z[mg.core_nodes] += U * dt

    # Test the results
    s = mg.at_node['topographic__steepest_slope']
    sa_factor = (1.0 + vs) * U / K
    a18 = mg.at_node['drainage_area'][18]
    a28 = mg.at_node['drainage_area'][28]
    s = mg.at_node['topographic__steepest_slope']
    s18 = sa_factor * (a18 ** -0.5)
    s28 = sa_factor * (a28 ** -0.5)
    assert_equal(np.round(s[18], 3), np.round(s18, 3))
    assert_equal(np.round(s[28], 3), np.round(s28, 3))

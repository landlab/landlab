# -*- coding: utf-8 -*-
"""
Unit tests for source_tracking_algorithm.py
@author: Sai Nudurupati & Erkan Istanbulluoglu
"""

from nose.tools import (assert_equal, assert_true, assert_raises, with_setup,
                        assert_almost_equal)
from numpy.testing import assert_array_almost_equal
try:
    from nose.tools import assert_is_instance
except ImportError:
    from landlab.testing.tools import assert_is_instance
import numpy as np

from landlab import RasterModelGrid
from landlab.components import FlowRouter
from landlab.utils import (track_source,
                           find_unique_upstream_hsd_ids_and_fractions)

def test_track_source():
    grid = RasterModelGrid((5, 5), spacing=(1., 1.))
    grid.at_node['topographic__elevation'] = np.array([5., 5., 5., 5., 5.,
                                                       5., 4., 5., 1., 5.,
                                                       0., 3., 5., 3., 0.,
                                                       5., 4., 5., 2., 5.,
                                                       5., 5., 5., 5., 5.])
    grid.status_at_node[10] = 0
    grid.status_at_node[14] = 0
    fr = FlowRouter(grid)
    fr.route_flow()
    r = grid.at_node['flow__receiver_node']
    assert_equal(r[6], 10)
    assert_equal(r[7], 8)
    assert_equal(r[18], 14)
    hsd_ids = np.empty(grid.number_of_nodes, dtype=int)
    hsd_ids[:] = 1
    hsd_ids[2:5] = 0
    hsd_ids[7:10] = 0
    (hsd_upstr, flow_accum) = track_source(grid, hsd_ids)
    assert_equal(hsd_upstr[8], [1, 0, 0])
    assert_equal(hsd_upstr[14], [1, 1, 1, 1, 0, 0, 1])
    assert_equal(flow_accum[14], 7)

def test_find_unique_upstream_hsd_ids_and_fractions():
    grid = RasterModelGrid((5, 5), spacing=(1., 1.))
    grid.at_node['topographic__elevation'] = np.array([5., 5., 5., 5., 5.,
                                                       5., 4., 5., 1., 5.,
                                                       0., 3., 5., 3., 0.,
                                                       5., 4., 5., 2., 5.,
                                                       5., 5., 5., 5., 5.])
    grid.status_at_node[10] = 0
    grid.status_at_node[14] = 0
    fr = FlowRouter(grid)
    fr.route_flow()
    hsd_ids = np.empty(grid.number_of_nodes, dtype=int)
    hsd_ids[:] = 1
    hsd_ids[2:5] = 0
    hsd_ids[7:10] = 0
    (hsd_upstr, flow_accum) = track_source(grid, hsd_ids)
    (uniq_ids, coeff) = find_unique_upstream_hsd_ids_and_fractions(hsd_upstr)
    assert_almost_equal(coeff[8][0], 0.666666667)
    
    
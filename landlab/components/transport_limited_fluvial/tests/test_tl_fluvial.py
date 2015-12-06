"""Test the transport limited fluvial module.

Simple test to ensure transport limited fluvial module runs and gives same
answers it always has.
"""
import os

import numpy as np
from numpy.testing import assert_array_equal, assert_array_almost_equal

from landlab import RasterModelGrid
from landlab.components.flow_routing.route_flow_dn import FlowRouter
from landlab.components.transport_limited_fluvial.tl_fluvial_monodirectional \
    import TransportLimitedEroder
from landlab import ModelParameterDictionary


_THIS_DIR = os.path.abspath(os.path.dirname(__file__))


def test_tl_fluvial():
    input_file = os.path.join(_THIS_DIR, 'stream_power_params_ideal.txt')
    inputs = ModelParameterDictionary(input_file)
    nrows = inputs.read_int('nrows')
    ncols = inputs.read_int('ncols')
    dx = inputs.read_float('dx')
    leftmost_elev = inputs.read_float('leftmost_elevation')
    initial_slope = inputs.read_float('initial_slope')
    uplift_rate = inputs.read_float('uplift_rate')

    runtime = inputs.read_float('total_time')
    dt = inputs.read_float('dt')

    nt = int(runtime // dt)
    uplift_per_step = uplift_rate * dt

    mg = RasterModelGrid(nrows, ncols, dx)
    mg.add_zeros('topographic__elevation', at='node')
    z = np.loadtxt(os.path.join(_THIS_DIR, 'tl_init.gz'))
    mg['node']['topographic__elevation'] = z

    mg.set_closed_boundaries_at_grid_edges(False, True, False, True)
    mg.set_fixed_value_boundaries_at_grid_edges(True, False, True, False,
                                                value_of='topographic__elevation')

    fr = FlowRouter(mg)
    tl = TransportLimitedEroder(mg, input_file)

    for i in xrange(nt):
        mg.at_node['topographic__elevation'][mg.core_nodes] += uplift_per_step
        mg = fr.route_flow()
        mg, _ = tl.erode(mg, dt, stability_condition='loose')

    z_tg = np.loadtxt(os.path.join(_THIS_DIR, 'tlz_tg.gz'))
    assert_array_almost_equal(mg.at_node['topographic__elevation'], z_tg)

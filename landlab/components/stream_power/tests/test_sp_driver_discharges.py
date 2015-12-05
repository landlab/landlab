"""Test simple stream power functionality when a discharge array is specified."""
import os

from six.moves import range

import numpy
import numpy as np
from numpy.testing import assert_array_equal, assert_array_almost_equal

from landlab import RasterModelGrid
from landlab import ModelParameterDictionary
from landlab.components.flow_routing.route_flow_dn import FlowRouter
from landlab.components.stream_power.stream_power import StreamPowerEroder


_THIS_DIR = os.path.abspath(os.path.dirname(__file__))


def test_sp_discharges():
    input_str = os.path.join(_THIS_DIR, 'test_sp_params_discharge.txt')
    inputs = ModelParameterDictionary(input_str)
    nrows = 5
    ncols = 5
    dx = inputs.read_float('dx')
    dt = inputs.read_float('dt')

    mg = RasterModelGrid(nrows, ncols, dx)
    mg.create_node_array_zeros('topographic__elevation')
    z = np.array([5., 5., 0., 5., 5.,
                  5., 2., 1., 2., 5.,
                  5., 3., 2., 3., 5.,
                  5., 4., 4., 4., 5.,
                  5., 5., 5., 5., 5.])
    mg['node']['topographic__elevation'] = z

    fr = FlowRouter(mg)
    sp = StreamPowerEroder(mg, input_str)

    # perform the loop (once!)
    for i in range(1):
        fr.route_flow(method='D8')
        my_Q = mg.at_node['water__volume_flux'] * 1.
        sp.erode(mg, dt, node_drainage_areas='drainage_area',
                 slopes_at_nodes='topographic__steepest_slope',
                 Q_if_used=my_Q)

    z_tg = np.array([5.00000000e+00,   5.00000000e+00,   0.00000000e+00,
                     5.00000000e+00,   5.00000000e+00,   5.00000000e+00,
                     1.29289322e+00,   1.00000000e-06,   1.29289322e+00,
                     5.00000000e+00,   5.00000000e+00,   2.29289322e+00,
                     1.00000000e+00,   2.29289322e+00,   5.00000000e+00,
                     5.00000000e+00,   3.29289322e+00,   3.00000000e+00,
                     3.29289322e+00,   5.00000000e+00,   5.00000000e+00,
                     5.00000000e+00,   5.00000000e+00,   5.00000000e+00,
                     5.00000000e+00])

    assert_array_almost_equal(mg.at_node['topographic__elevation'], z_tg)

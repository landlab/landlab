"""Test simple stream power functionality when a width array is specified."""
import os

from six.moves import range

import numpy
import numpy as np
from numpy.testing import assert_array_equal, assert_array_almost_equal

from landlab import RasterModelGrid
from landlab import ModelParameterDictionary
from landlab.components.flow_routing import FlowRouter
from landlab.components.stream_power import StreamPowerEroder


_THIS_DIR = os.path.abspath(os.path.dirname(__file__))


def test_sp_widths():
    input_str = os.path.join(_THIS_DIR, 'test_sp_params_widths.txt')
    inputs = ModelParameterDictionary(input_str, auto_type=True)
    nrows = 5
    ncols = 5
    dx = inputs.read_float('dx')
    dt = inputs.read_float('dt')

    mg = RasterModelGrid(nrows, ncols, dx)
    widths = np.ones(mg.number_of_nodes, dtype=float)
    mg.add_zeros('topographic__elevation', at='node')
    z = np.array([5., 5., 0., 5., 5.,
                  5., 2., 1., 2., 5.,
                  5., 3., 2., 3., 5.,
                  5., 4., 4., 4., 5.,
                  5., 5., 5., 5., 5.])
    mg['node']['topographic__elevation'] = z

    fr = FlowRouter(mg)
    sp = StreamPowerEroder(mg, use_W=widths, **inputs)

    # perform the loop (once!)
    for i in range(1):
        fr.route_flow()
        sqrt_A = mg.at_node['drainage_area']**0.5
        widths[mg.core_nodes] = sqrt_A[mg.core_nodes]/sqrt_A[
            mg.core_nodes].mean()
        # so widths has mean=1.
        # note the issue with drainage_area not defined at perimeter => nans
        # if not careful...
        sp.run_one_step(dt)

    z_tg = np.array([ 5.        ,  5.        ,  0.        ,  5.        ,
                      5.        ,  5.        ,  1.37222369,  0.36876358,
                      1.37222369,  5.        ,  5.        ,  2.17408606,
                      1.07986038,  2.17408606,  5.        ,  5.        ,
                      3.08340277,  2.85288049,  3.08340277,  5.        ,
                      5.        ,  5.        ,  5.        ,  5.        ,
                      5.        ])

    assert_array_almost_equal(mg.at_node['topographic__elevation'], z_tg)
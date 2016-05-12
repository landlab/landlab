"""Test the storm generator and simple stream power eroder both execute."""
import os
import pylab
import time
import copy

import numpy as np
from numpy.testing import assert_array_equal, assert_array_almost_equal

try:
    from nose.tools import assert_is
except ImportError:
    from landlab.testing.tools import assert_is

from landlab import RasterModelGrid
from landlab import ModelParameterDictionary
from landlab.components.flow_routing import FlowRouter
from landlab.components.stream_power import StreamPowerEroder
from landlab.components.uniform_precip import PrecipitationDistribution


_THIS_DIR = os.path.abspath(os.path.dirname(__file__))


def test_storms():
    input_file_string = os.path.join(_THIS_DIR, 'drive_sp_params_storms.txt')
    inputs = ModelParameterDictionary(input_file_string, auto_type=True)
    nrows = inputs.read_int('nrows')
    ncols = inputs.read_int('ncols')
    dx = inputs.read_float('dx')
    dt = inputs.read_float('dt')
    time_to_run = inputs.read_float('run_time')
    uplift = inputs.read_float('uplift_rate')

    mean_duration  = inputs.read_float('mean_storm')
    mean_interstorm  = inputs.read_float('mean_interstorm')
    mean_depth = inputs.read_float('mean_depth')

    storm_run_time  = inputs.read_float('storm_run_time')
    delta_t  = inputs.read_float('delta_t')
    mg = RasterModelGrid(nrows, ncols, dx)

    mg.add_zeros('topographic__elevation', at='node')
    z = mg.zeros(at='node')
    mg['node']['topographic__elevation'] = z + np.random.rand(len(z)) / 1000.
    mg.add_zeros('water__unit_flux_in', at='node')

    precip = PrecipitationDistribution(mean_storm_duration = mean_duration,
                                mean_interstorm_duration = mean_interstorm,
                                mean_storm_depth = mean_depth,
                                total_t = storm_run_time, delta_t = delta_t)
    fr = FlowRouter(mg)
    sp = StreamPowerEroder(mg, **inputs)

    for (interval_duration, rainfall_rate) in \
            precip.yield_storm_interstorm_duration_intensity():
        if rainfall_rate != 0.:
            mg.at_node['water__unit_flux_in'].fill(rainfall_rate)
            mg = fr.route_flow()
            sp.run_one_step(dt)
        mg.at_node['topographic__elevation'][
            mg.core_nodes] += uplift * interval_duration

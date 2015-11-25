'''
Test the storm generator and simple stream power eroder both execute.
'''
from numpy.testing import assert_array_equal, assert_array_almost_equal
try:
    from nose.tools import assert_is
except ImportError:
    from landlab.testing.tools import assert_is
from landlab.components.flow_routing.route_flow_dn import FlowRouter
from landlab.components.stream_power.stream_power import StreamPowerEroder
from landlab.components.uniform_precip.generate_uniform_precip import \
    PrecipitationDistribution

import numpy
from landlab import RasterModelGrid
from landlab import ModelParameterDictionary
import pylab
import time
import copy


def test_storms():
    input_file_string = ('../landlab/components/stream_power/tests/' +
                         'drive_sp_params_storms.txt')
    inputs = ModelParameterDictionary(input_file_string)
    nrows = inputs.read_int('nrows')
    ncols = inputs.read_int('ncols')
    dx = inputs.read_float('dx')
    dt = inputs.read_float('dt')
    time_to_run = inputs.read_float('run_time')
    uplift = inputs.read_float('uplift_rate')

    mg = RasterModelGrid(nrows, ncols, dx)

    mg.create_node_array_zeros('topographic__elevation')
    z = mg.create_node_array_zeros()
    mg['node']['topographic__elevation'] = z + numpy.random.rand(len(z))/1000.
    mg.add_zeros('node', 'water__volume_flux_in')

    precip = PrecipitationDistribution(input_file=input_file_string)
    fr = FlowRouter(mg)
    sp = StreamPowerEroder(mg, input_file_string)

    for (interval_duration, rainfall_rate) in \
            precip.yield_storm_interstorm_duration_intensity():
        if rainfall_rate != 0.:
            mg.at_node['water__volume_flux_in'].fill(rainfall_rate)
            mg = fr.route_flow()
            sp.erode(mg, dt=interval_duration, Q_if_used='water__volume_flux')
        mg.at_node['topographic__elevation'][
            mg.core_nodes] += uplift * interval_duration

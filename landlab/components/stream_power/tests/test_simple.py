"""Test Braun-Willett and non-fastscape stream power component.

A simple driver implementing Braun-Willett flow routing and then a
(non-fastscape) stream power component.
"""
# DEJH, 09/15/14

import os

import numpy
from numpy.testing import assert_array_almost_equal

from landlab import RasterModelGrid
from landlab import ModelParameterDictionary
from landlab.components.flow_routing.route_flow_dn import FlowRouter
from landlab.components.stream_power.stream_power import StreamPowerEroder


_THIS_DIR = os.path.abspath(os.path.dirname(__file__))


def test_fastscape():
    input_str = os.path.join(_THIS_DIR, 'drive_sp_params.txt')
    inputs = ModelParameterDictionary(input_str)
    nrows = inputs.read_int('nrows')
    ncols = inputs.read_int('ncols')
    dx = inputs.read_float('dx')
    dt = inputs.read_float('dt')
    time_to_run = inputs.read_float('run_time')
    uplift = inputs.read_float('uplift_rate')
    init_elev = inputs.read_float('init_elev')

    mg = RasterModelGrid((nrows, ncols), spacing=(dx, dx))
    mg.set_closed_boundaries_at_grid_edges(False, False, True, True)

    mg.create_node_array_zeros('topographic__elevation')
    z = mg.create_node_array_zeros() + init_elev
    numpy.random.seed(0)
    mg['node']['topographic__elevation'] = z + \
        numpy.random.rand(len(z)) / 1000.

    fr = FlowRouter(mg)
    sp = StreamPowerEroder(mg, input_str)
    elapsed_time = 0.
    while elapsed_time < time_to_run:
        if elapsed_time + dt > time_to_run:
            dt = time_to_run - elapsed_time
        fr.route_flow(method='D8')
        sp.erode(mg, dt)
        mg.at_node['topographic__elevation'][mg.core_nodes] += uplift * dt
        elapsed_time += dt

    z_trg = numpy.array([5.48813504e-04,   7.15189366e-04,   6.02763376e-04,
                         5.44883183e-04,   4.23654799e-04,   6.45894113e-04,
                         1.02783376e-02,   9.66667235e-03,   6.15060782e-03,
                         3.83441519e-04,   7.91725038e-04,   1.00905776e-02,
                         8.98955843e-03,   5.32836181e-03,   7.10360582e-05,
                         8.71292997e-05,   9.96377080e-03,   9.63738797e-03,
                         7.31213677e-03,   8.70012148e-04,   9.78618342e-04,
                         1.02124693e-02,   8.78386002e-03,   4.88161060e-03,
                         1.18274426e-04,   6.39921021e-04,   1.00377580e-02,
                         9.54340293e-03,   6.05173814e-03,   4.14661940e-04,
                         2.64555612e-04,   1.02160196e-02,   8.61600088e-03,
                         4.77005225e-03,   1.87898004e-05,   6.17635497e-04,
                         9.30445558e-03,   9.48713993e-03,   7.25689742e-03,
                         6.81820299e-04,   3.59507901e-04,   5.79161813e-03,
                         6.83777542e-03,   6.18063842e-03,   6.66766715e-04,
                         6.70637870e-04,   2.10382561e-04,   1.28926298e-04,
                         3.15428351e-04,   3.63710771e-04])

    assert_array_almost_equal(mg.at_node['topographic__elevation'], z_trg)

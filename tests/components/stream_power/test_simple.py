"""Test Braun-Willett and non-fastscape stream power component.

A simple driver implementing Braun-Willett flow routing and then a
(non-fastscape) stream power component.
"""

# DEJH, 09/15/14

import os

import numpy
from numpy.testing import assert_array_almost_equal

from landlab import RasterModelGrid
from landlab.components import FlowAccumulator
from landlab.components import StreamPowerEroder

_THIS_DIR = os.path.abspath(os.path.dirname(__file__))


def test_sp_old():
    uplift = 0.001
    dt = 1.0
    time_to_run = 10.0
    mg = RasterModelGrid((10, 5))
    mg.set_closed_boundaries_at_grid_edges(False, False, True, True)

    mg.add_zeros("topographic__elevation", at="node")
    z = mg.zeros(at="node")
    numpy.random.seed(0)
    mg["node"]["topographic__elevation"] = z + numpy.random.rand(len(z)) / 1000.0

    fr = FlowAccumulator(mg, flow_director="D8")
    sp = StreamPowerEroder(mg, K_sp=0.1, m_sp=0.5, n_sp=1.0, threshold_sp=0.0)
    elapsed_time = 0.0
    while elapsed_time < time_to_run:
        if elapsed_time + dt > time_to_run:
            dt = time_to_run - elapsed_time
        fr.run_one_step()
        sp.run_one_step(dt)
        mg.at_node["topographic__elevation"][mg.core_nodes] += uplift * dt
        elapsed_time += dt

    z_trg = numpy.array(
        [
            5.48813504e-04,
            7.15189366e-04,
            6.02763376e-04,
            5.44883183e-04,
            4.23654799e-04,
            6.45894113e-04,
            1.01830760e-02,
            9.58036770e-03,
            6.55865452e-03,
            3.83441519e-04,
            7.91725038e-04,
            1.00142749e-02,
            8.80798884e-03,
            5.78387585e-03,
            7.10360582e-05,
            8.71292997e-05,
            9.81911417e-03,
            9.52243406e-03,
            7.55093226e-03,
            8.70012148e-04,
            9.78618342e-04,
            1.00629755e-02,
            8.49253798e-03,
            5.33216680e-03,
            1.18274426e-04,
            6.39921021e-04,
            9.88956320e-03,
            9.47119567e-03,
            6.43790696e-03,
            4.14661940e-04,
            2.64555612e-04,
            1.00450743e-02,
            8.37262908e-03,
            5.21540904e-03,
            1.87898004e-05,
            6.17635497e-04,
            9.21286940e-03,
            9.34022513e-03,
            7.51114450e-03,
            6.81820299e-04,
            3.59507901e-04,
            6.19166921e-03,
            7.10456176e-03,
            6.62585507e-03,
            6.66766715e-04,
            6.70637870e-04,
            2.10382561e-04,
            1.28926298e-04,
            3.15428351e-04,
            3.63710771e-04,
        ]
    )

    assert_array_almost_equal(mg.at_node["topographic__elevation"], z_trg)


def test_sp_new():
    """
    Tests new style component instantiation and run.
    """
    uplift = 0.001
    dt = 1.0
    time_to_run = 10.0
    mg = RasterModelGrid((10, 5))

    mg.set_closed_boundaries_at_grid_edges(False, False, True, True)

    mg.add_zeros("topographic__elevation", at="node")
    z = mg.zeros(at="node")
    numpy.random.seed(0)
    mg["node"]["topographic__elevation"] = z + numpy.random.rand(len(z)) / 1000.0

    fr = FlowAccumulator(mg, flow_director="D8")
    sp = StreamPowerEroder(mg, K_sp=0.1, m_sp=0.5, n_sp=1.0, threshold_sp=0.0)
    elapsed_time = 0.0
    while elapsed_time < time_to_run:
        if elapsed_time + dt > time_to_run:
            dt = time_to_run - elapsed_time
        fr.run_one_step()
        sp.run_one_step(dt)
        mg.at_node["topographic__elevation"][mg.core_nodes] += uplift * dt
        elapsed_time += dt

    z_trg = numpy.array(
        [
            5.48813504e-04,
            7.15189366e-04,
            6.02763376e-04,
            5.44883183e-04,
            4.23654799e-04,
            6.45894113e-04,
            1.01830760e-02,
            9.58036770e-03,
            6.55865452e-03,
            3.83441519e-04,
            7.91725038e-04,
            1.00142749e-02,
            8.80798884e-03,
            5.78387585e-03,
            7.10360582e-05,
            8.71292997e-05,
            9.81911417e-03,
            9.52243406e-03,
            7.55093226e-03,
            8.70012148e-04,
            9.78618342e-04,
            1.00629755e-02,
            8.49253798e-03,
            5.33216680e-03,
            1.18274426e-04,
            6.39921021e-04,
            9.88956320e-03,
            9.47119567e-03,
            6.43790696e-03,
            4.14661940e-04,
            2.64555612e-04,
            1.00450743e-02,
            8.37262908e-03,
            5.21540904e-03,
            1.87898004e-05,
            6.17635497e-04,
            9.21286940e-03,
            9.34022513e-03,
            7.51114450e-03,
            6.81820299e-04,
            3.59507901e-04,
            6.19166921e-03,
            7.10456176e-03,
            6.62585507e-03,
            6.66766715e-04,
            6.70637870e-04,
            2.10382561e-04,
            1.28926298e-04,
            3.15428351e-04,
            3.63710771e-04,
        ]
    )

    assert_array_almost_equal(mg.at_node["topographic__elevation"], z_trg)

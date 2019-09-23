"""Test Braun-Willett flow routing and fastscape stream power.

A simple driver implementing Braun-Willett flow routing and then the
fastscape stream power component.
"""
# DEJH, 09/15/14
import os

import numpy
from numpy.testing import assert_array_almost_equal

from landlab import ModelParameterDictionary, RasterModelGrid
from landlab.components import FlowAccumulator
from landlab.components.stream_power import FastscapeEroder as Fsc

_THIS_DIR = os.path.abspath(os.path.dirname(__file__))


def test_fastscape():
    input_str = os.path.join(_THIS_DIR, "drive_sp_params.txt")
    inputs = ModelParameterDictionary(input_str)
    nrows = inputs.read_int("nrows")
    ncols = inputs.read_int("ncols")
    dx = inputs.read_float("dx")
    dt = inputs.read_float("dt")
    time_to_run = inputs.read_float("run_time")
    uplift = inputs.read_float("uplift_rate")
    init_elev = inputs.read_float("init_elev")

    mg = RasterModelGrid((nrows, ncols), xy_spacing=dx)
    mg.set_closed_boundaries_at_grid_edges(False, False, True, True)

    mg.add_zeros("topographic__elevation", at="node")
    z = mg.zeros(at="node") + init_elev
    numpy.random.seed(0)
    mg["node"]["topographic__elevation"] = z + numpy.random.rand(len(z)) / 1000.0

    fr = FlowAccumulator(mg, flow_director="D8")
    fsp = Fsc(mg, input_str, method="D8")
    elapsed_time = 0.0
    while elapsed_time < time_to_run:
        if elapsed_time + dt > time_to_run:
            dt = time_to_run - elapsed_time
        mg = fr.run_one_step()
        mg = fsp.erode(mg, dt=dt)
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


def test_fastscape_new():
    input_str = os.path.join(_THIS_DIR, "drive_sp_params.txt")
    inputs = ModelParameterDictionary(input_str, auto_type=True)
    nrows = inputs.read_int("nrows")
    ncols = inputs.read_int("ncols")
    dx = inputs.read_float("dx")
    dt = inputs.read_float("dt")
    time_to_run = inputs.read_float("run_time")
    uplift = inputs.read_float("uplift_rate")
    init_elev = inputs.read_float("init_elev")

    mg = RasterModelGrid((nrows, ncols), xy_spacing=dx)
    mg.set_closed_boundaries_at_grid_edges(False, False, True, True)

    mg.add_zeros("topographic__elevation", at="node")
    z = mg.zeros(at="node") + init_elev
    numpy.random.seed(0)
    mg["node"]["topographic__elevation"] = z + numpy.random.rand(len(z)) / 1000.0

    fr = FlowAccumulator(mg, flow_director="D8")
    fsp = Fsc(mg, **inputs)  # here's the diff from the above
    elapsed_time = 0.0
    while elapsed_time < time_to_run:
        if elapsed_time + dt > time_to_run:
            dt = time_to_run - elapsed_time
        fr.run_one_step()
        fsp.run_one_step(dt)  # new style
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

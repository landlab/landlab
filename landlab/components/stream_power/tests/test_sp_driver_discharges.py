"""Test simple stream power functionality when a discharge array is specified."""
import os

import numpy as np
from numpy.testing import assert_array_almost_equal
from six.moves import range

from landlab import ModelParameterDictionary, RasterModelGrid
from landlab.components import FlowAccumulator, StreamPowerEroder

_THIS_DIR = os.path.abspath(os.path.dirname(__file__))


def test_sp_discharges_old():
    input_str = os.path.join(_THIS_DIR, "test_sp_params_discharge.txt")
    inputs = ModelParameterDictionary(input_str)
    nrows = 5
    ncols = 5
    dx = inputs.read_float("dx")
    dt = inputs.read_float("dt")

    mg = RasterModelGrid(nrows, ncols, dx)
    mg.add_zeros("topographic__elevation", at="node")
    z = np.array(
        [
            5.,
            5.,
            0.,
            5.,
            5.,
            5.,
            2.,
            1.,
            2.,
            5.,
            5.,
            3.,
            2.,
            3.,
            5.,
            5.,
            4.,
            4.,
            4.,
            5.,
            5.,
            5.,
            5.,
            5.,
            5.,
        ]
    )
    mg["node"]["topographic__elevation"] = z

    fr = FlowAccumulator(mg, flow_director="D8")
    my_Q = mg.at_node["surface_water__discharge"]
    sp = StreamPowerEroder(mg, input_str, use_Q=my_Q)

    # perform the loop (once!)
    for i in range(1):
        fr.run_one_step()
        my_Q[:] = mg.at_node["surface_water__discharge"] * 1.
        sp.run_one_step(dt)

    z_tg = np.array(
        [
            5.,
            5.,
            0.,
            5.,
            5.,
            5.,
            1.47759225,
            0.43050087,
            1.47759225,
            5.,
            5.,
            2.32883687,
            1.21525044,
            2.32883687,
            5.,
            5.,
            3.27261262,
            3.07175015,
            3.27261262,
            5.,
            5.,
            5.,
            5.,
            5.,
            5.,
        ]
    )

    assert_array_almost_equal(mg.at_node["topographic__elevation"], z_tg)


def test_sp_discharges_new():
    input_str = os.path.join(_THIS_DIR, "test_sp_params_discharge_new.txt")
    inputs = ModelParameterDictionary(input_str, auto_type=True)
    nrows = 5
    ncols = 5
    dx = inputs.read_float("dx")
    dt = inputs.read_float("dt")

    mg = RasterModelGrid(nrows, ncols, dx)
    mg.add_zeros("topographic__elevation", at="node")
    z = np.array(
        [
            5.,
            5.,
            0.,
            5.,
            5.,
            5.,
            2.,
            1.,
            2.,
            5.,
            5.,
            3.,
            2.,
            3.,
            5.,
            5.,
            4.,
            4.,
            4.,
            5.,
            5.,
            5.,
            5.,
            5.,
            5.,
        ]
    )
    mg["node"]["topographic__elevation"] = z

    fr = FlowAccumulator(mg, flow_director="D8")
    sp = StreamPowerEroder(mg, **inputs)

    # perform the loop (once!)
    for i in range(1):
        fr.run_one_step()
        sp.run_one_step(dt)

    z_tg = np.array(
        [
            5.,
            5.,
            0.,
            5.,
            5.,
            5.,
            1.47759225,
            0.43050087,
            1.47759225,
            5.,
            5.,
            2.32883687,
            1.21525044,
            2.32883687,
            5.,
            5.,
            3.27261262,
            3.07175015,
            3.27261262,
            5.,
            5.,
            5.,
            5.,
            5.,
            5.,
        ]
    )

    assert_array_almost_equal(mg.at_node["topographic__elevation"], z_tg)

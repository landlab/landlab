"""Test simple stream power functionality when a width array is specified."""
import os

import numpy as np
from numpy.testing import assert_array_almost_equal
from six.moves import range

from landlab import ModelParameterDictionary, RasterModelGrid
from landlab.components import FlowAccumulator, StreamPowerEroder

_THIS_DIR = os.path.abspath(os.path.dirname(__file__))


def test_sp_widths():
    input_str = os.path.join(_THIS_DIR, "test_sp_params_widths.txt")
    inputs = ModelParameterDictionary(input_str, auto_type=True)
    nrows = 5
    ncols = 5
    dx = inputs.read_float("dx")
    dt = inputs.read_float("dt")

    mg = RasterModelGrid((nrows, ncols), xy_spacing=dx)
    widths = np.ones(mg.number_of_nodes, dtype=float)
    mg.add_zeros("topographic__elevation", at="node")
    z = np.array(
        [
            5.0,
            5.0,
            0.0,
            5.0,
            5.0,
            5.0,
            2.0,
            1.0,
            2.0,
            5.0,
            5.0,
            3.0,
            2.0,
            3.0,
            5.0,
            5.0,
            4.0,
            4.0,
            4.0,
            5.0,
            5.0,
            5.0,
            5.0,
            5.0,
            5.0,
        ]
    )
    mg["node"]["topographic__elevation"] = z

    fr = FlowAccumulator(mg, flow_director="D8")
    sp = StreamPowerEroder(mg, use_W=widths, **inputs)

    # perform the loop (once!)
    for i in range(1):
        fr.run_one_step()
        sqrt_A = mg.at_node["drainage_area"] ** 0.5
        widths[mg.core_nodes] = sqrt_A[mg.core_nodes] / sqrt_A[mg.core_nodes].mean()
        # so widths has mean=1.
        # note the issue with drainage_area not defined at perimeter => nans
        # if not careful...
        sp.run_one_step(dt)

    z_tg = np.array(
        [
            5.0,
            5.0,
            0.0,
            5.0,
            5.0,
            5.0,
            1.37222369,
            0.36876358,
            1.37222369,
            5.0,
            5.0,
            2.17408606,
            1.07986038,
            2.17408606,
            5.0,
            5.0,
            3.08340277,
            2.85288049,
            3.08340277,
            5.0,
            5.0,
            5.0,
            5.0,
            5.0,
            5.0,
        ]
    )

    assert_array_almost_equal(mg.at_node["topographic__elevation"], z_tg)

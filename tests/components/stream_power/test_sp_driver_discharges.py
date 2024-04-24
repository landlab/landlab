"""Test simple stream power functionality when a discharge array is specified."""

import os

import numpy as np
from numpy.testing import assert_array_almost_equal

from landlab import RasterModelGrid
from landlab.components import FlowAccumulator
from landlab.components import StreamPowerEroder

_THIS_DIR = os.path.abspath(os.path.dirname(__file__))


def test_sp_discharges_old():
    dt = 1.0
    mg = RasterModelGrid((5, 5))
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
    sp = StreamPowerEroder(
        mg, K_sp=0.5, m_sp=0.5, n_sp=1.0, discharge_field="surface_water__discharge"
    )

    # perform the loop (once!)
    for _ in range(1):
        fr.run_one_step()
        sp.run_one_step(dt)

    z_tg = np.array(
        [
            5.0,
            5.0,
            0.0,
            5.0,
            5.0,
            5.0,
            1.47759225,
            0.43050087,
            1.47759225,
            5.0,
            5.0,
            2.32883687,
            1.21525044,
            2.32883687,
            5.0,
            5.0,
            3.27261262,
            3.07175015,
            3.27261262,
            5.0,
            5.0,
            5.0,
            5.0,
            5.0,
            5.0,
        ]
    )

    assert_array_almost_equal(mg.at_node["topographic__elevation"], z_tg)


def test_sp_discharges_new():
    dt = 1.0

    mg = RasterModelGrid((5, 5))
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
    sp = StreamPowerEroder(mg, K_sp=0.5, m_sp=0.5, n_sp=1.0, threshold_sp=0.0)

    # perform the loop (once!)
    for _ in range(1):
        fr.run_one_step()
        sp.run_one_step(dt)

    z_tg = np.array(
        [
            5.0,
            5.0,
            0.0,
            5.0,
            5.0,
            5.0,
            1.47759225,
            0.43050087,
            1.47759225,
            5.0,
            5.0,
            2.32883687,
            1.21525044,
            2.32883687,
            5.0,
            5.0,
            3.27261262,
            3.07175015,
            3.27261262,
            5.0,
            5.0,
            5.0,
            5.0,
            5.0,
            5.0,
        ]
    )

    assert_array_almost_equal(mg.at_node["topographic__elevation"], z_tg)

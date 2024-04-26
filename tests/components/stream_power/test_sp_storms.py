"""Test the storm generator and simple stream power eroder both execute."""

import os

import numpy as np

from landlab import RasterModelGrid
from landlab.components import FlowAccumulator
from landlab.components import PrecipitationDistribution
from landlab.components import StreamPowerEroder

_THIS_DIR = os.path.abspath(os.path.dirname(__file__))


def test_storms():
    dt = 500.0
    uplift = 0.0001

    mean_duration = 100.0
    mean_interstorm = 400.0
    mean_depth = 5.0

    storm_run_time = 3000000.0
    delta_t = 500.0
    mg = RasterModelGrid((10, 10), xy_spacing=1000.0)

    mg.add_zeros("topographic__elevation", at="node")
    z = mg.zeros(at="node")
    mg["node"]["topographic__elevation"] = z + np.random.rand(len(z)) / 1000.0
    mg.add_zeros("water__unit_flux_in", at="node")

    precip = PrecipitationDistribution(
        mean_storm_duration=mean_duration,
        mean_interstorm_duration=mean_interstorm,
        mean_storm_depth=mean_depth,
        total_t=storm_run_time,
        delta_t=delta_t,
    )
    fr = FlowAccumulator(mg, flow_director="D8")
    sp = StreamPowerEroder(
        mg,
        K_sp=0.0001,
        m_sp=0.5,
        n_sp=1.0,
        threshold_sp=1.0,
        discharge_field="surface_water__discharge",
    )

    for (
        interval_duration,
        rainfall_rate,
    ) in precip.yield_storm_interstorm_duration_intensity():
        if rainfall_rate != 0.0:
            mg.at_node["water__unit_flux_in"].fill(rainfall_rate)
            fr.run_one_step()
            sp.run_one_step(dt)
        mg.at_node["topographic__elevation"][mg.core_nodes] += (
            uplift * interval_duration
        )

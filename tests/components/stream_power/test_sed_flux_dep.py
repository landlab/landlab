"""Test the SedDepEroder component.

Test the sed dep eroder by turning it over a few times. No attempt has been
made to ensure the solution is stable. Takes a topo already output and runs it
a few more times, to ensure repeatability.
"""
import os

import numpy as np
from numpy.testing import assert_array_almost_equal

from landlab import RasterModelGrid
from landlab.components import FlowAccumulator, SedDepEroder

_THIS_DIR = os.path.abspath(os.path.dirname(__file__))


def test_sed_dep():
    uplift_rate = 0.0001
    runtime = 20000.0
    dt = 5000.0

    nt = int(runtime // dt)
    uplift_per_step = uplift_rate * dt

    mg = RasterModelGrid((50, 50), xy_spacing=1000.0)

    mg.add_zeros("topographic__elevation", at="node")
    z = np.loadtxt(os.path.join(_THIS_DIR, "seddepinit.txt"))
    mg["node"]["topographic__elevation"] = z

    mg.set_closed_boundaries_at_grid_edges(True, False, True, False)

    fr = FlowAccumulator(mg, flow_director="D8")
    sde = SedDepEroder(
        mg,
        K_sp=1.0e-6,
        b_sp=0.5,
        c_sp=1.0,
        Qc="MPM",
        k_Q=2.5e-7,
        k_w=2.5,
        mannings_n=0.05,
        Dchar=0.05,
        threshold_Shields=0.05,
        sed_dependency_type="almost_parabolic",
        set_threshold_from_Dchar=True,
        g=9.81,
    )

    for i in range(nt):
        mg.at_node["topographic__elevation"][mg.core_nodes] += uplift_per_step
        mg = fr.run_one_step()
        mg, _ = sde.run_one_step(dt)

    z_tg = np.loadtxt(os.path.join(_THIS_DIR, "seddepz_tg.txt"))

    assert_array_almost_equal(
        mg.at_node["topographic__elevation"][mg.core_nodes], z_tg[mg.core_nodes]
    )


def test_sed_dep_new():
    """
    This tests only the power_law version of the SDE.
    It uses a landscape run to 5000 100y iterations, then having experienced
    a 20-fold uplift acceleration for a further 30000 y. It tests the outcome
    of the next 1000 y of erosion.
    """
    mg = RasterModelGrid((25, 50), xy_spacing=200.0)
    for edge in (mg.nodes_at_left_edge, mg.nodes_at_top_edge, mg.nodes_at_right_edge):
        mg.status_at_node[edge] = mg.BC_NODE_IS_CLOSED

    z = mg.add_zeros("topographic__elevation", at="node")

    fr = FlowAccumulator(mg, flow_director="D8")
    sde = SedDepEroder(
        mg,
        K_sp=1.0e-4,
        sed_dependency_type="almost_parabolic",
        Qc="power_law",
        K_t=1.0e-4,
    )

    initconds = os.path.join(os.path.dirname(__file__), "perturbedcondst300.txt")
    finalconds = os.path.join(os.path.dirname(__file__), "tenmorestepsfrom300.txt")
    z[:] = np.loadtxt(initconds)

    dt = 100.0
    up = 0.05

    for i in range(10):
        fr.run_one_step()
        sde.run_one_step(dt)
        z[mg.core_nodes] += 20.0 * up

    assert_array_almost_equal(z, np.loadtxt(finalconds))

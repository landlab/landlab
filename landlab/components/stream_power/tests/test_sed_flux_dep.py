"""Test the SedDepEroder component.

Test the sed dep eroder by turning it over a few times. No attempt has been
made to ensure the solution is stable. Takes a topo already output and runs it
a few more times, to ensure repeatability.
"""
import os

from six.moves import range

import numpy as np
import os
from numpy.testing import assert_array_almost_equal

from landlab import RasterModelGrid, CLOSED_BOUNDARY
from landlab.components.flow_routing import FlowRouter
from landlab.components.stream_power import SedDepEroder
from landlab import ModelParameterDictionary


_THIS_DIR = os.path.abspath(os.path.dirname(__file__))


def test_sed_dep_new():
    """
    This tests only the power_law version of the SDE.
    It uses a landscape run to 5000 100y iterations, then having experienced
    a 20-fold uplift acceleration for a further 30000 y. It tests the outcome
    of the next 1000 y of erosion.
    """
    mg = RasterModelGrid((25, 50), 200.)
    for edge in (mg.nodes_at_left_edge, mg.nodes_at_top_edge,
                 mg.nodes_at_right_edge):
        mg.status_at_node[edge] = CLOSED_BOUNDARY

    z = mg.add_zeros('node', 'topographic__elevation')

    fr = FlowRouter(mg)
    sde = SedDepEroder(mg, K_sp=1.e-4, sed_dependency_type='almost_parabolic',
                       Qc='power_law', K_t=1.e-4)

    initconds = os.path.join(os.path.dirname(__file__),
                             'perturbedcondst300.txt')
    finalconds = os.path.join(os.path.dirname(__file__),
                              'tenmorestepsfrom300.txt')
    z[:] = np.loadtxt(initconds)

    dt = 100.
    up = 0.05

    for i in range(10):
        fr.route_flow()
        sde.run_one_step(dt)
        z[mg.core_nodes] += 20.*up

    assert_array_almost_equal(z, np.loadtxt(finalconds))

"""Test the SedDepEroder component.

Test the sed dep eroder by turning it over a few times. No attempt has been
made to ensure the solution is stable. Takes a topo already output and runs it
a few more times, to ensure repeatability.
"""
import os

from six.moves import range

import numpy as np
import os
from numpy.testing import assert_array_almost_equal, assert_equal

from landlab import RasterModelGrid, CLOSED_BOUNDARY
from landlab.components.flow_routing import FlowRouter
from landlab.components.stream_power import SedDepEroder
from landlab import ModelParameterDictionary


def test_sed_dep_new():
    """
    This tests only the power_law version of the SDE.
    """
    mg = RasterModelGrid((10, 5), 200.)
    for edge in (mg.nodes_at_left_edge, mg.nodes_at_top_edge,
                 mg.nodes_at_right_edge):
        mg.status_at_node[edge] = CLOSED_BOUNDARY

    z = mg.add_zeros('node', 'topographic__elevation')

    fr = FlowRouter(mg)
    sde = SedDepEroder(mg, K_sp=1.e-4, sed_dependency_type='almost_parabolic',
                       Qc='power_law', K_t=1.e-4, external_sediment=True)

    z[:] = mg.node_y/10000.
    z.reshape((10, 5))[:, 2] *= 2.
    z.reshape((10, 5))[:, 3] *= 2.1

    initz = z.copy()

    dt = 100.
    up = 0.05

    for i in range(1):
        print(i)
        fr.route_flow()
        sde.run_one_step(dt)

    assert_array_almost_equal(mg.at_node['channel_sediment__depth'],
                              np.array([0.,  0.,  0.,  0.,  0.,
                                        0.,  0.,  0.,  0.,  0.,
                                        0.,  0.,  0.,  0.,  0.,
                                        0.,  0.,  0.,  0.,  0.,
                                        0.,  0.,  0.,  0.,  0.,
                                        0.,  0.,  0.,  0.,  0.,
                                        0.,  0.,  0.,  0.,  0.,
                                        0.,  0.,  0.,  0.,  0.,
                                        0.,  0.00023431,  0.,  0.,  0.,
                                        0.,  0.,  0.,  0.,  0.]))

    assert_array_almost_equal(z - initz, np.array(
        [0.,  0.,  0.,  0.,  0.,
         0., -0.00076662, -0.0004,     -0.00118507,  0.,
         0., -0.00068314, -0.00042426, -0.00110741,  0.,
         0., -0.00065571, -0.0006,     -0.00102346,  0.,
         0., -0.00056135, -0.0008,     -0.00093111,  0.,
         0., -0.00045180, -0.0010,     -0.00078882,  0.,
         0., -0.00031903, -0.0012,     -0.00071743,  0.,
         0., -0.00015763, -0.0014,     -0.00057541,  0.,
         0.,  0.00023431, -0.0016,     -0.00042,     0.,
         0.,  0.,  0.,  0.,  0.]))

    assert_equal(len(sde._error_at_abort), 0)  # good convergence at all nodes


def test_sed_dep_w_hillslopes():
    """
    This tests only the power_law version of the SDE, with a hillslope input.
    """
    mg = RasterModelGrid((10, 5), 200.)
    for edge in (mg.nodes_at_left_edge, mg.nodes_at_top_edge,
                 mg.nodes_at_right_edge):
        mg.status_at_node[edge] = CLOSED_BOUNDARY

    z = mg.add_zeros('node', 'topographic__elevation')
    th = mg.add_zeros('node', 'channel_sediment__depth')
    th[mg.core_nodes] += 0.001

    fr = FlowRouter(mg)
    sde = SedDepEroder(mg, K_sp=1.e-4, sed_dependency_type='almost_parabolic',
                       Qc='power_law', K_t=1.e-4, external_sediment=True)

    z[:] = mg.node_y/10000.
    z.reshape((10, 5))[:, 2] *= 2.
    z.reshape((10, 5))[:, 3] *= 2.1

    initz = z.copy()

    dt = 100.
    up = 0.05

    for i in range(1):
        print(i)
        fr.route_flow()
        sde.run_one_step(dt)

    # test binding of field occurs correctly:
    assert th is mg.at_node['channel_sediment__depth']

    assert_array_almost_equal(
        mg.at_node['channel_sediment__depth'],
        np.array([0.00000000e+00,   0.00000000e+00,   0.00000000e+00,
                  0.00000000e+00,   0.00000000e+00,   0.00000000e+00,
                  0.00000000e+00,   6.00000000e-04,   0.00000000e+00,
                  0.00000000e+00,   0.00000000e+00,   0.00000000e+00,
                  5.75735931e-04,   0.00000000e+00,   0.00000000e+00,
                  0.00000000e+00,   0.00000000e+00,   4.00000000e-04,
                  0.00000000e+00,   0.00000000e+00,   0.00000000e+00,
                  9.28079257e-07,   2.00000000e-04,   0.00000000e+00,
                  0.00000000e+00,   0.00000000e+00,   4.13904292e-04,
                  0.00000000e+00,   0.00000000e+00,   0.00000000e+00,
                  0.00000000e+00,   7.60612309e-04,   0.00000000e+00,
                  5.55537486e-06,   0.00000000e+00,   0.00000000e+00,
                  1.16568542e-03,   0.00000000e+00,   2.32060608e-04,
                  0.00000000e+00,   0.00000000e+00,   1.73431458e-03,
                  0.00000000e+00,   5.80000000e-04,   0.00000000e+00,
                  0.00000000e+00,   0.00000000e+00,   0.00000000e+00,
                  0.00000000e+00,   0.00000000e+00]))

    assert_array_almost_equal(z, np.array(
        [0.,  0.,  0.,  0.,  0.,
         0.02,  0.01877957,  0.0396,      0.04055596,  0.02,
         0.04,  0.03891440,  0.07957574,  0.08263689,  0.04,
         0.06,  0.05890177,  0.1194,      0.12472345,  0.06,
         0.08,  0.07900093,  0.1592,      0.16681590,  0.08,
         0.10,  0.09941390,  0.1990,      0.20891231,  0.10,
         0.12,  0.11976061,  0.23863333,  0.25100556,  0.12,
         0.14,  0.14016569,  0.27831429,  0.29323206,  0.14,
         0.16,  0.16073431,  0.318025,    0.335580,    0.16,
         0.18,  0.18,  0.36,  0.378,  0.18]))

    # good convergence at all nodes
    assert_equal(len(sde._error_at_abort), 0)

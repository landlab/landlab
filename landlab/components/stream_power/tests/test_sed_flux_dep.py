"""Test the SedDepEroder component.

Test the sed dep eroder by turning it over a few times. No attempt has been
made to ensure the solution is stable. Takes a topo already output and runs it
a few more times, to ensure repeatability.
"""
import os

import numpy as np
import os
from numpy.testing import assert_array_almost_equal, assert_equal
from six.moves import range

from landlab import RasterModelGrid, CLOSED_BOUNDARY, ModelParameterDictionary
from landlab.components import FlowAccumulator
from landlab.components import SedDepEroder
from landlab.components import FastscapeEroder


def test_sed_dep_new_almostpara():
    """
    This tests only the power_law version of the SDE, using the
    almost_parabolic form of f(Qs).
    """
    mg = RasterModelGrid((10, 5), xy_spacing=200.)
    for edge in (mg.nodes_at_left_edge, mg.nodes_at_top_edge,
                 mg.nodes_at_right_edge):
        mg.status_at_node[edge] = CLOSED_BOUNDARY

    z = mg.add_zeros('node', 'topographic__elevation')

    fr = FlowAccumulator(mg, flow_director='D8')
    sde = SedDepEroder(mg, K_sp=1.e-4, sed_dependency_type='almost_parabolic',
                       Qc='power_law', K_t=1.e-4)

    z[:] = mg.node_y/10000.
    z.reshape((10, 5))[:, 2] *= 2.
    z.reshape((10, 5))[:, 3] *= 2.1

    initz = z.copy()

    dt = 100.
    up = 0.05

    for i in range(1):
        fr.run_one_step()
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


def test_sed_dep():
    """
    This may be now invalid, and has got accidentally reintroduced during
    conflict resolution.
    """
    input_file = os.path.join(_THIS_DIR, "sed_dep_params.txt")
    inputs = ModelParameterDictionary(input_file, auto_type=True)
    nrows = inputs.read_int("nrows")
    ncols = inputs.read_int("ncols")
    dx = inputs.read_float("dx")
    uplift_rate = inputs.read_float("uplift_rate")

    runtime = inputs.read_float("total_time")
    dt = inputs.read_float("dt")

    initz = z.copy()

    mg = RasterModelGrid((nrows, ncols), xy_spacing=(dx, dx))

    mg.add_zeros("topographic__elevation", at="node")
    z = np.loadtxt(os.path.join(_THIS_DIR, "seddepinit.txt"))
    mg["node"]["topographic__elevation"] = z


def test_sed_dep_new_genhumped():
    """
    This tests only the power_law version of the SDE, using the
    generalized_humped form of f(Qs).
    """
    mg = RasterModelGrid((10, 5), xy_spacing=200.)
    for edge in (mg.nodes_at_left_edge, mg.nodes_at_top_edge,
                 mg.nodes_at_right_edge):
        mg.status_at_node[edge] = CLOSED_BOUNDARY

    z = mg.add_zeros('node', 'topographic__elevation')

    fr = FlowAccumulator(mg, flow_director='D8')
    sde = SedDepEroder(mg, K_sp=1.e-4,
                       sed_dependency_type='generalized_humped',
                       Qc='power_law', K_t=1.e-4)

    z[:] = mg.node_y/10000.
    z.reshape((10, 5))[:, 2] *= 2.
    z.reshape((10, 5))[:, 3] *= 2.1

    initz = z.copy()

    dt = 100.
    up = 0.05

    for i in range(1):
        fr.run_one_step()
        sde.run_one_step(dt)

    ans = np.array([0.,  0.,  0.,  0.,  0.,
                    0.02,  0.01940896,  0.03969863,  0.04112046,  0.02,
                    0.04,  0.03948890,  0.07968035,  0.08318039,  0.04,
                    0.06,  0.05951858,  0.11954794,  0.12524486,  0.06,
                    0.08,  0.07960747,  0.15939726,  0.16731501,  0.08,
                    0.10,  0.09969946,  0.19924657,  0.20939250,  0.10,
                    0.12,  0.11979225,  0.23909589,  0.25147973,  0.12,
                    0.14,  0.13987971,  0.27894520,  0.29357949,  0.14,
                    0.16,  0.16023431,  0.31879451,  0.33568356,  0.16,
                    0.18,  0.18,  0.36,  0.378,  0.18])

    assert_array_almost_equal(z, ans)


def test_sed_dep_new_lindecl():
    """
    This tests only the power_law version of the SDE, using the
    linear_decline form of f(Qs).
    """
    mg = RasterModelGrid((10, 5), xy_spacing=200.)
    for edge in (
        mg.nodes_at_left_edge, mg.nodes_at_top_edge, mg.nodes_at_right_edge
    ):
        mg.status_at_node[edge] = CLOSED_BOUNDARY

    z = mg.add_zeros("node", "topographic__elevation")

    fr = FlowAccumulator(mg, flow_director='D8')
    sde = SedDepEroder(
        mg,
        K_sp=1.e-4,
        sed_dependency_type='linear_decline',
        Qc='power_law',
        K_t=1.e-4)

    z[:] = mg.node_y/10000.
    z.reshape((10, 5))[:, 2] *= 2.
    z.reshape((10, 5))[:, 3] *= 2.1

    initz = z.copy()

    dt = 100.
    up = 0.05

    for i in range(1):
        fr.run_one_step()
        sde.run_one_step(dt)

    ans = np.array([0.,  0.,  0.,  0.,  0.,
                    0.02,  0.01955879,  0.03980000,  0.04128996,  0.02,
                    0.04,  0.03961955,  0.07978787,  0.08333878,  0.04,
                    0.06,  0.05964633,  0.11970000,  0.12539158,  0.06,
                    0.08,  0.07971138,  0.15960000,  0.16745207,  0.08,
                    0.10,  0.09978034,  0.19950000,  0.20951474,  0.10,
                    0.12,  0.11985392,  0.23940000,  0.25158663,  0.12,
                    0.14,  0.13993079,  0.27930000,  0.29367338,  0.14,
                    0.16,  0.16023431,  0.31920000,  0.33579000,  0.16,
                    0.18,  0.18,  0.36,  0.378,  0.18])

    assert_array_almost_equal(z, ans)


def test_sed_dep_new_const():
    """
    This tests only the power_law version of the SDE, using the
    constant (None) form of f(Qs).
    """
    mg = RasterModelGrid((10, 5), xy_spacing=200.)
    for edge in (
        mg.nodes_at_left_edge, mg.nodes_at_top_edge, mg.nodes_at_right_edge
    ):
        mg.status_at_node[edge] = CLOSED_BOUNDARY

    z = mg.add_zeros('node', 'topographic__elevation')

    fr = FlowAccumulator(mg, flow_director='D8')
    sde = SedDepEroder(
        mg,
        K_sp=1.e-4,
        sed_dependency_type='None',
        Qc='power_law',
        K_t=1.e-4
    )

    z[:] = mg.node_y/10000.
    z.reshape((10, 5))[:, 2] *= 2.
    z.reshape((10, 5))[:, 3] *= 2.1

    initz = z.copy()

    dt = 100.
    up = 0.05

    for i in range(1):
        fr.run_one_step()
        sde.run_one_step(dt)

    ans = np.array([0.,  0.,  0.,  0.,  0.,
                    0.02,  0.01922540,  0.03960000,  0.04081206,  0.02,
                    0.04,  0.03927889,  0.07957574,  0.08288878,  0.04,
                    0.06,  0.05930718,  0.11940000,  0.12497121,  0.06,
                    0.08,  0.07936754,  0.15920000,  0.16706085,  0.08,
                    0.10,  0.09943431,  0.19900000,  0.20916000,  0.10,
                    0.12,  0.11951010,  0.23880000,  0.25127254,  0.12,
                    0.14,  0.13960000,  0.27860000,  0.29340603,  0.14,
                    0.16,  0.16023431,  0.31840000,  0.33558000,  0.16,
                    0.18,  0.18,  0.36,  0.378,  0.18])

    assert_array_almost_equal(z, ans)


def test_sed_dep_w_hillslopes():
    """
    This tests only the power_law version of the SDE, with a hillslope input.
    """
    mg = RasterModelGrid((10, 5), xy_spacing=200.)
    for edge in (
        mg.nodes_at_left_edge, mg.nodes_at_top_edge, mg.nodes_at_right_edge
    ):
        mg.status_at_node[edge] = CLOSED_BOUNDARY

    z = mg.add_zeros('node', 'topographic__elevation')
    th = mg.add_zeros('node', 'channel_sediment__depth')
    th[mg.core_nodes] += 0.001

    fr = FlowAccumulator(mg, flow_director='D8')
    sde = SedDepEroder(
        mg,
        K_sp=1.e-4,
        sed_dependency_type='almost_parabolic',
        Qc='power_law',
        K_t=1.e-4
    )

    z[:] = mg.node_y/10000.
    z.reshape((10, 5))[:, 2] *= 2.
    z.reshape((10, 5))[:, 3] *= 2.1

    initz = z.copy()

    dt = 100.
    up = 0.05

    for i in range(1):
        print(i)
        fr.run_one_step()
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

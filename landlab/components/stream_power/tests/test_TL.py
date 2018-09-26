import pytest
import numpy as np

from numpy.testing import assert_array_almost_equal

from landlab import RasterModelGrid
from landlab.components import TransportLimitedEroder
from landlab.components import FlowAccumulator


def test_multidirection():
    mg = RasterModelGrid(4, 4, 1)
    z = mg.add_zeros('node', 'topographic__elevation')
    fa = FlowAccumulator(mg, flow_director='D8')
    fa.run_one_step()
    with pytest.raises(NotImplementedError):
        TL = TransportLimitedEroder(mg)


def test_bad_init1():
    mg = RasterModelGrid(4, 4, 1)
    with pytest.raises(ValueError):
        TL = TransportLimitedEroder(mg, phi=1.)


def test_bad_init2():
    mg = RasterModelGrid(4, 4, 1)
    with pytest.raises(ValueError):
        TL = TransportLimitedEroder(mg, phi=-1.)


def test_bad_init3():
    mg = RasterModelGrid(4, 4, 1)
    with pytest.raises(ValueError):
        TL = TransportLimitedEroder(mg, F_f=1.1)


def test_bad_init4():
    mg = RasterModelGrid(4, 4, 1)
    with pytest.raises(ValueError):
        TL = TransportLimitedEroder(mg, F_f=-1.)


def test_bad_init5():
    mg = RasterModelGrid(4, 4, 1)
    with pytest.raises(ValueError):
        TL = TransportLimitedEroder(mg, solver='nope')


def test_bad_init6():
    mg = RasterModelGrid(4, 4, 1)
    with pytest.raises(ValueError):
        TL = TransportLimitedEroder(mg, sp_crit=-1.)


def test_bad_init7():
    mg = RasterModelGrid(4, 4, 1)
    badcrits = np.zeros(16, dtype=float)
    badcrits[5] = -1.
    with pytest.raises(ValueError):
        TL = TransportLimitedEroder(mg, sp_crit=badcrits)


def test_bad_init8():
    mg = RasterModelGrid(4, 4, 1)
    badcrits = np.zeros(16, dtype=float)
    badcrits[5] = -1.
    crits = mg.add_field('node', 'crits', badcrits)
    with pytest.raises(ValueError):
        TL = TransportLimitedEroder(mg, sp_crit='crits')


def test_sed_flux_in():
    mg = RasterModelGrid(4, 4, 1)
    Qsin = mg.add_ones('node', 'sediment__flux')
    TL = TransportLimitedEroder(mg)
    assert np.allclose(TL.qs, 1.)






# def test_no_thresh():
#     K = 0.001
#     U = 0.01
#     m = 0.5
#     n = 1.0
#     threshold = 0.0
#     dt = 1000
#
#     mg = RasterModelGrid(30, 3, 100.)
#     mg.set_closed_boundaries_at_grid_edges(True, False, True, False)
#     z = mg.zeros(at='node')
#     mg['node']['topographic__elevation'] = z + np.random.rand(len(z)) / 1000.
#
#     fa = FlowAccumulator(mg)
#     sp = Spst(mg, K_sp = K, threshold_sp=threshold)
#     for i in range(100):
#         fa.run_one_step()
#         sp.run_one_step(dt)
#         mg['node']['topographic__elevation'][mg.core_nodes] += U*dt
#
#     actual_slopes = mg.at_node["topographic__steepest_slope"][mg.core_nodes[1:-1]]
#     actual_areas = mg.at_node["drainage_area"][mg.core_nodes[1:-1]]
#
#     predicted_slopes = (U / (K * (actual_areas ** m))) ** (1. / n)
#
#     assert_array_almost_equal(actual_slopes, predicted_slopes)
#
#
# def test_with_thresh():
#     K = 0.001
#     U = 0.01
#     m = 0.5
#     n = 1.0
#     threshold = 1.0
#     dt = 1000
#
#     mg = RasterModelGrid(30, 3, 100.)
#     mg.set_closed_boundaries_at_grid_edges(True, False, True, False)
#     z = mg.zeros(at='node')
#     mg['node']['topographic__elevation'] = z + np.random.rand(len(z)) / 1000.
#
#     fa = FlowAccumulator(mg)
#     sp = Spst(mg, K_sp = K, threshold_sp=threshold)
#     for i in range(100):
#         fa.run_one_step()
#         sp.run_one_step(dt)
#         mg['node']['topographic__elevation'][mg.core_nodes] += U*dt
#
#     actual_slopes = mg.at_node["topographic__steepest_slope"][mg.core_nodes[1:-1]]
#     actual_areas = mg.at_node["drainage_area"][mg.core_nodes[1:-1]]
#
#     predicted_slopes_upper = ((U + threshold) / (K * (actual_areas ** m))) ** (1. / n)
#     predicted_slopes_lower = ((U + 0.0) / (K * (actual_areas ** m))) ** (1. / n)
#
#     # assert actual and predicted slopes are in the correct range for the slopes.
#     assert np.all(actual_slopes > predicted_slopes_lower) == True
#     assert np.all(actual_slopes < predicted_slopes_upper) == True

import numpy as np
import pytest
from numpy.testing import assert_array_almost_equal

from landlab import RasterModelGrid
from landlab.components import FlowAccumulator, StreamPowerSmoothThresholdEroder as Spst


def test_bad_nsp():

    mg = RasterModelGrid((4, 4))
    mg.set_closed_boundaries_at_grid_edges(False, False, True, True)
    with pytest.raises(ValueError):
        Spst(mg, K_sp=1.0, n_sp=1.01)


def test_no_thresh():
    K = 0.001
    U = 0.01
    m = 0.5
    n = 1.0
    threshold = 0.0
    dt = 1000

    mg = RasterModelGrid((30, 3), xy_spacing=100.0)
    mg.set_closed_boundaries_at_grid_edges(True, False, True, False)
    z = mg.zeros(at="node")
    mg["node"]["topographic__elevation"] = z + np.random.rand(len(z)) / 1000.0

    fa = FlowAccumulator(mg)
    sp = Spst(mg, K_sp=K, threshold_sp=threshold)
    for i in range(100):
        fa.run_one_step()
        sp.run_one_step(dt)
        mg["node"]["topographic__elevation"][mg.core_nodes] += U * dt

    actual_slopes = mg.at_node["topographic__steepest_slope"][mg.core_nodes[1:-1]]
    actual_areas = mg.at_node["drainage_area"][mg.core_nodes[1:-1]]

    predicted_slopes = (U / (K * (actual_areas ** m))) ** (1.0 / n)

    assert_array_almost_equal(actual_slopes, predicted_slopes)


def test_with_thresh():
    K = 0.001
    U = 0.01
    m = 0.5
    n = 1.0
    threshold = 1.0
    dt = 1000

    mg = RasterModelGrid((30, 3), xy_spacing=100.0)
    mg.set_closed_boundaries_at_grid_edges(True, False, True, False)
    z = mg.zeros(at="node")
    mg["node"]["topographic__elevation"] = z + np.random.rand(len(z)) / 1000.0

    fa = FlowAccumulator(mg)
    sp = Spst(mg, K_sp=K, threshold_sp=threshold)
    for i in range(100):
        fa.run_one_step()
        sp.run_one_step(dt)
        mg["node"]["topographic__elevation"][mg.core_nodes] += U * dt

    actual_slopes = mg.at_node["topographic__steepest_slope"][mg.core_nodes[1:-1]]
    actual_areas = mg.at_node["drainage_area"][mg.core_nodes[1:-1]]

    predicted_slopes_upper = ((U + threshold) / (K * (actual_areas ** m))) ** (1.0 / n)
    predicted_slopes_lower = ((U + 0.0) / (K * (actual_areas ** m))) ** (1.0 / n)

    # assert actual and predicted slopes are in the correct range for the slopes.
    assert np.all(actual_slopes > predicted_slopes_lower)
    assert np.all(actual_slopes < predicted_slopes_upper)

import numpy as np
import pytest
from numpy.testing import assert_array_equal

from landlab import RasterModelGrid
from landlab.components import DrainageDensity, FastscapeEroder, FlowAccumulator


def test_route_to_multiple_error_raised():
    mg = RasterModelGrid((10, 10))
    z = mg.add_zeros("node", "topographic__elevation")
    z += mg.x_of_node + mg.y_of_node
    fa = FlowAccumulator(mg, flow_director="MFD")
    fa.run_one_step()

    channel__mask = mg.zeros(at="node")

    with pytest.raises(NotImplementedError):
        DrainageDensity(mg, channel__mask=channel__mask)


def test_mask_is_stable():
    mg = RasterModelGrid((80, 80), 1.0)
    mg.add_zeros("node", "topographic__elevation")
    np.random.seed(50)
    noise = np.random.rand(mg.size("node"))
    mg.at_node["topographic__elevation"] += noise
    mg.at_node["topographic__elevation"]  # doctest: +NORMALIZE_WHITESPACE
    fr = FlowAccumulator(mg, flow_director="D8")
    fsc = FastscapeEroder(mg, K_sp=.01, m_sp=.5, n_sp=1)
    for x in range(100):
        fr.run_one_step()
        fsc.run_one_step(dt=10.0)
        mg.at_node["topographic__elevation"][mg.core_nodes] += .01

    mask = np.zeros(len(mg.at_node["topographic__elevation"]), dtype=np.uint8)
    mask[np.where(mg.at_node["drainage_area"] > 5)] = 1

    mask0 = mask.copy()

    dd = DrainageDensity(mg, channel__mask=mask)
    mask1 = mask.copy()

    dd.calc_drainage_density()
    mask2 = mask.copy()

    assert_array_equal(mask0, mask1)
    assert_array_equal(mask0[mg.core_nodes], mask2[mg.core_nodes])

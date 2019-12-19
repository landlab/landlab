import numpy as np
import pytest
from numpy.testing import assert_array_equal

from landlab import FieldError, RasterModelGrid
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
    mg = RasterModelGrid((10, 10))
    mg.add_zeros("node", "topographic__elevation")
    np.random.seed(3542)
    noise = np.random.rand(mg.size("node"))
    mg.at_node["topographic__elevation"] += noise
    fr = FlowAccumulator(mg, flow_director="D8")
    fsc = FastscapeEroder(mg, K_sp=0.01, m_sp=0.5, n_sp=1)
    for x in range(2):
        fr.run_one_step()
        fsc.run_one_step(dt=10.0)
        mg.at_node["topographic__elevation"][mg.core_nodes] += 0.01

    mask = np.zeros(len(mg.at_node["topographic__elevation"]), dtype=np.uint8)
    mask[np.where(mg.at_node["drainage_area"] > 0)] = 1

    mask0 = mask.copy()

    dd = DrainageDensity(mg, channel__mask=mask)
    mask1 = mask.copy()

    dd.calc_drainage_density()
    mask2 = mask.copy()

    assert_array_equal(mask0, mask1)
    assert_array_equal(mask0[mg.core_nodes], mask2[mg.core_nodes])


def test_float_mask():
    mg = RasterModelGrid((10, 10))
    mg.add_zeros("node", "topographic__elevation")
    noise = np.random.rand(mg.size("node"))
    mg.at_node["topographic__elevation"] += noise
    fr = FlowAccumulator(mg, flow_director="D8")
    fr.run_one_step()
    mask = np.zeros(len(mg.at_node["topographic__elevation"]), dtype=float)
    mask[np.where(mg.at_node["drainage_area"] > 5)] = 1
    with pytest.raises(ValueError):
        DrainageDensity(mg, channel__mask=mask)


def test_bool_mask():
    mg = RasterModelGrid((10, 10))
    mg.add_zeros("node", "topographic__elevation")
    noise = np.random.rand(mg.size("node"))
    mg.at_node["topographic__elevation"] += noise
    fr = FlowAccumulator(mg, flow_director="D8")
    fr.run_one_step()
    mask = np.zeros(len(mg.at_node["topographic__elevation"]), dtype=bool)
    mask[np.where(mg.at_node["drainage_area"] > 5)] = 1
    with pytest.raises(ValueError):
        DrainageDensity(mg, channel__mask=mask)


def test_missing_fields():
    mg = RasterModelGrid((10, 10))
    mg.add_zeros("node", "topographic__elevation")
    with pytest.raises(FieldError):
        DrainageDensity(
            mg,
            area_coefficient=1,
            slope_coefficient=1,
            area_exponent=1,
            slope_exponent=1,
            channelization_threshold=1,
        )


def test_updating_with_array_provided():
    mg = RasterModelGrid((10, 10))
    mg.add_zeros("node", "topographic__elevation")
    noise = np.random.rand(mg.size("node"))
    mg.at_node["topographic__elevation"] += noise
    fr = FlowAccumulator(mg, flow_director="D8")
    fr.run_one_step()
    mask = np.zeros(len(mg.at_node["topographic__elevation"]), dtype=np.uint8)
    mask[np.where(mg.at_node["drainage_area"] > 5)] = 1
    dd = DrainageDensity(mg, channel__mask=mask)
    with pytest.raises(NotImplementedError):
        dd._update_channel_mask()


def test_mask_field_exists():
    mg = RasterModelGrid((10, 10))
    mg.add_zeros("node", "topographic__elevation")
    mg.add_zeros("node", "channel__mask")
    noise = np.random.rand(mg.size("node"))
    mg.at_node["topographic__elevation"] += noise
    fr = FlowAccumulator(mg, flow_director="D8")
    fr.run_one_step()
    mask = np.zeros(len(mg.at_node["topographic__elevation"]), dtype=np.uint8)
    mask[np.where(mg.at_node["drainage_area"] > 5)] = 1
    with pytest.warns(UserWarning):
        DrainageDensity(mg, channel__mask=mask)


def test_bad_mask_size():
    mg = RasterModelGrid((10, 10))
    mg.add_zeros("node", "topographic__elevation")
    noise = np.random.rand(mg.size("node"))
    mg.at_node["topographic__elevation"] += noise
    fr = FlowAccumulator(mg, flow_director="D8")
    fr.run_one_step()
    mask = np.zeros(20)
    mask[np.where(mg.at_node["drainage_area"][:20] > 5)] = 1

    with pytest.raises(ValueError):
        DrainageDensity(mg, channel__mask=mask)


def test_providing_array_and_kwargs():
    mg = RasterModelGrid((10, 10))
    mg.add_zeros("node", "topographic__elevation")
    noise = np.random.rand(mg.size("node"))
    mg.at_node["topographic__elevation"] += noise
    fr = FlowAccumulator(mg, flow_director="D8")
    fr.run_one_step()
    mask = np.zeros(len(mg.at_node["topographic__elevation"]), dtype=np.uint8)
    mask[np.where(mg.at_node["drainage_area"] > 5)] = 1
    with pytest.warns(UserWarning):
        DrainageDensity(mg, channel__mask=mask, area_coefficient=1)
    with pytest.warns(UserWarning):
        DrainageDensity(mg, channel__mask=mask, slope_coefficient=1)
    with pytest.warns(UserWarning):
        DrainageDensity(mg, channel__mask=mask, area_exponent=1)
    with pytest.warns(UserWarning):
        DrainageDensity(mg, channel__mask=mask, slope_exponent=1)
    with pytest.warns(UserWarning):
        DrainageDensity(mg, channel__mask=mask, channelization_threshold=1)


def test_missing_channel_threshold():
    mg = RasterModelGrid((10, 10))
    mg.add_zeros("node", "topographic__elevation")
    noise = np.random.rand(mg.size("node"))
    mg.at_node["topographic__elevation"] += noise
    fr = FlowAccumulator(mg, flow_director="D8")
    fr.run_one_step()
    with pytest.raises(ValueError):
        DrainageDensity(
            mg,
            area_coefficient=1,
            slope_coefficient=1,
            area_exponent=1,
            slope_exponent=1,
        )


def test_missing_slope_exponent():
    mg = RasterModelGrid((10, 10))
    mg.add_zeros("node", "topographic__elevation")
    noise = np.random.rand(mg.size("node"))
    mg.at_node["topographic__elevation"] += noise
    fr = FlowAccumulator(mg, flow_director="D8")
    fr.run_one_step()
    with pytest.raises(ValueError):
        DrainageDensity(
            mg,
            area_coefficient=1,
            slope_coefficient=1,
            area_exponent=1,
            channelization_threshold=1,
        )


def test_missing_slope_coefficient():
    mg = RasterModelGrid((10, 10))
    mg.add_zeros("node", "topographic__elevation")
    noise = np.random.rand(mg.size("node"))
    mg.at_node["topographic__elevation"] += noise
    fr = FlowAccumulator(mg, flow_director="D8")
    fr.run_one_step()
    with pytest.raises(ValueError):
        DrainageDensity(
            mg,
            area_coefficient=1,
            area_exponent=1,
            slope_exponent=1,
            channelization_threshold=1,
        )


def test_missing_area_exponent():
    mg = RasterModelGrid((10, 10))
    mg.add_zeros("node", "topographic__elevation")
    noise = np.random.rand(mg.size("node"))
    mg.at_node["topographic__elevation"] += noise
    fr = FlowAccumulator(mg, flow_director="D8")
    fr.run_one_step()
    with pytest.raises(ValueError):
        DrainageDensity(
            mg,
            area_coefficient=1,
            slope_coefficient=1,
            slope_exponent=1,
            channelization_threshold=1,
        )


def test_missing_area_coefficient():
    mg = RasterModelGrid((10, 10))
    mg.add_zeros("node", "topographic__elevation")
    noise = np.random.rand(mg.size("node"))
    mg.at_node["topographic__elevation"] += noise
    fr = FlowAccumulator(mg, flow_director="D8")
    fr.run_one_step()
    with pytest.raises(ValueError):
        DrainageDensity(
            mg,
            slope_coefficient=1,
            area_exponent=1,
            slope_exponent=1,
            channelization_threshold=1,
        )

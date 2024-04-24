import pytest

from landlab import RasterModelGrid
from landlab.components import FastscapeEroder
from landlab.components import FlowAccumulator
from landlab.components import SedDepEroder
from landlab.components import StreamPowerEroder
from landlab.components import StreamPowerSmoothThresholdEroder


def test_route_to_multiple_error_raised_init_FastscapeEroder():
    mg = RasterModelGrid((10, 10))
    z = mg.add_zeros("node", "topographic__elevation")
    z += mg.x_of_node + mg.y_of_node
    fa = FlowAccumulator(mg, flow_director="MFD")
    fa.run_one_step()

    with pytest.raises(NotImplementedError):
        FastscapeEroder(mg)


def test_route_to_multiple_error_raised_init_SedDepEroder():
    mg = RasterModelGrid((10, 10))
    z = mg.add_zeros("node", "topographic__elevation")
    z += mg.x_of_node + mg.y_of_node
    fa = FlowAccumulator(mg, flow_director="MFD")
    fa.run_one_step()

    with pytest.raises(NotImplementedError):
        SedDepEroder(mg)


def test_route_to_multiple_error_raised_init_StreamPowerSmoothThresholdEroder():
    mg = RasterModelGrid((10, 10))
    z = mg.add_zeros("node", "topographic__elevation")
    z += mg.x_of_node + mg.y_of_node
    fa = FlowAccumulator(mg, flow_director="MFD")
    fa.run_one_step()

    with pytest.raises(NotImplementedError):
        StreamPowerSmoothThresholdEroder(mg)


def test_route_to_multiple_error_raised_init_StreamPowerEroder():
    mg = RasterModelGrid((10, 10))
    z = mg.add_zeros("node", "topographic__elevation")
    z += mg.x_of_node + mg.y_of_node
    fa = FlowAccumulator(mg, flow_director="MFD")
    fa.run_one_step()

    with pytest.raises(NotImplementedError):
        StreamPowerEroder(mg)

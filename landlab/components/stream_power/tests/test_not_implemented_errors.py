import pytest
from landlab import RasterModelGrid
from landlab.components import (FlowAccumulator,
                                FastscapeEroder,
                                SedDepEroder,
                                StreamPowerSmoothThresholdEroder)


def test_route_to_multiple_error_raised_FastscapeEroder():
    mg = RasterModelGrid((10, 10))
    z = mg.add_zeros('node', 'topographic__elevation')
    z += mg.x_of_node + mg.y_of_node
    fa = FlowAccumulator(mg, flow_director='MFD')
    fa.run_one_step()

    with pytest.raises(NotImplementedError):
        FastscapeEroder(mg)


def test_route_to_multiple_error_raised_SedDepEroder():
    mg = RasterModelGrid((10, 10))
    z = mg.add_zeros('node', 'topographic__elevation')
    z += mg.x_of_node + mg.y_of_node
    fa = FlowAccumulator(mg, flow_director='MFD')
    fa.run_one_step()

    with pytest.raises(NotImplementedError):
        SedDepEroder(mg)


def test_route_to_multiple_error_raised_StreamPowerSmoothThresholdEroder():
    mg = RasterModelGrid((10, 10))
    z = mg.add_zeros('node', 'topographic__elevation')
    z += mg.x_of_node + mg.y_of_node
    fa = FlowAccumulator(mg, flow_director='MFD')
    fa.run_one_step()

    with pytest.raises(NotImplementedError):
        StreamPowerSmoothThresholdEroder(mg)


def test_route_to_multiple_error_raised_StreamPowerEroder():
    mg = RasterModelGrid((10, 10))
    z = mg.add_zeros('node', 'topographic__elevation')
    z += mg.x_of_node + mg.y_of_node
    fa = FlowAccumulator(mg, flow_director='MFD')
    fa.run_one_step()

    with pytest.raises(NotImplementedError):
        StreamPowerEroder(mg)

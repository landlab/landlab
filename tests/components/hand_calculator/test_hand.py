import numpy as np
import pytest

from landlab import RasterModelGrid
from landlab.components import FlowAccumulator
from landlab.components import HeightAboveDrainageCalculator


def test_route_to_multiple_error_raised():
    mg = RasterModelGrid((10, 10))
    z = mg.add_zeros("topographic__elevation", at="node")
    z += mg.x_of_node + mg.y_of_node
    fa = FlowAccumulator(mg, flow_director="MFD")
    fa.run_one_step()

    channel__mask = mg.zeros(at="node")

    with pytest.raises(NotImplementedError):
        HeightAboveDrainageCalculator(mg, channel_mask=channel__mask)


def test_warn_drainage_pits():
    mg = RasterModelGrid((4, 4))
    z = mg.add_zeros("topographic__elevation", at="node")
    elev = np.array([[2, 1, 1, 2], [3, 2, 2, 3], [4, 1, 3, 4], [5, 3, 4, 4]])
    z[:] = elev.reshape(len(z))

    fa = FlowAccumulator(mg, flow_director="D8")
    fa.run_one_step()

    channel__mask = mg.zeros(at="node")
    channel__mask[[2, 6]] = 1
    hd = HeightAboveDrainageCalculator(mg, channel_mask=channel__mask)

    with pytest.warns(UserWarning):
        hd.run_one_step()

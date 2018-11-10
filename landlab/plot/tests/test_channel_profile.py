import pytest

from landlab import RasterModelGrid
from landlab.components import FlowAccumulator
from landlab.plot import analyze_channel_network_and_plot
from landlab.plot.channel_profile import channel_nodes


def test_route_to_multiple_error_raised():
    mg = RasterModelGrid((10, 10))
    z = mg.add_zeros("node", "topographic__elevation")
    z += mg.x_of_node + mg.y_of_node
    fa = FlowAccumulator(mg, flow_director="MFD")
    fa.run_one_step()

    with pytest.raises(NotImplementedError):
        analyze_channel_network_and_plot(mg)

    with pytest.raises(NotImplementedError):
        channel_nodes(
            mg,
            mg.at_node["topographic__steepest_slope"],
            mg.at_node["drainage_area"],
            mg.at_node["flow__receiver_node"],
            number_of_channels=1,
            threshold=1.0,
        )

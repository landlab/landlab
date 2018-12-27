import pytest
from landlab import RasterModelGrid, NetworkModelGrid


@pytest.fixture
def four_by_four_raster():
    mg = RasterModelGrid((4, 4))
    return mg


@pytest.fixture
def simple_network():
    y_of_node = (0, 1, 2, 2)
    x_of_node = (0, 0, -1, 1)
    nodes_at_link = ((1, 0), (2, 1), (3, 1))
    mg = NetworkModelGrid((y_of_node, x_of_node), nodes_at_link)
    return mg

import pytest

from landlab import NetworkModelGrid, RasterModelGrid
from landlab.values.synthetic import _STATUS


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


def pytest_generate_tests(metafunc):
    if "at" in metafunc.fixturenames:
        metafunc.parametrize("at", ("node", "link", "patch", "corner", "face", "cell"))
    if "node_bc" in metafunc.fixturenames:
        metafunc.parametrize("node_bc", list(_STATUS["node"].keys()))
    if "link_bc" in metafunc.fixturenames:
        metafunc.parametrize("link_bc", list(_STATUS["link"].keys()))

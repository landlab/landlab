from landlab.graph.hex.hex import HorizontalHexTriGraph
from landlab.graph.hex.hex import HorizontalRectTriGraph
from landlab.graph.hex.hex import VerticalHexTriGraph
from landlab.graph.hex.hex import VerticalRectTriGraph


def pytest_generate_tests(metafunc):
    if "hex_layout" in metafunc.fixturenames:
        metafunc.parametrize(
            "hex_layout",
            [
                HorizontalRectTriGraph,
                VerticalRectTriGraph,
                HorizontalHexTriGraph,
                VerticalHexTriGraph,
            ],
        )

from landlab.graph.hex.hex import (
    HorizontalHexTriGraph,
    HorizontalRectTriGraph,
    VerticalHexTriGraph,
    VerticalRectTriGraph,
)


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

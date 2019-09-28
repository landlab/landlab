import numpy as np
import pytest

from landlab.graph.hex.hex import (
    HorizontalHexTriGraph,
    HorizontalRectTriGraph,
    VerticalHexTriGraph,
    VerticalRectTriGraph,
)


@pytest.fixture
def hex_shape():
    return tuple(np.random.randint(0, 1000, 2))


@pytest.fixture
def small_hex_shape():
    return tuple(np.random.randint(0, 100, 2))


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

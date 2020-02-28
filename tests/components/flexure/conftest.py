import pytest

from landlab import RasterModelGrid
from landlab.components.flexure.flexure import Flexure
from landlab.components.flexure.flexure_1d import Flexure1D


@pytest.fixture
def flex():
    grid = RasterModelGrid((20, 20), xy_spacing=10e3)
    grid.add_zeros("lithosphere__overlying_pressure_increment", at="node")
    return Flexure(grid)


@pytest.fixture
def flex1d():
    grid = RasterModelGrid((20, 20), xy_spacing=10e3)
    grid.add_zeros("lithosphere_surface__increment_of_elevation", at="node")
    grid.add_zeros("lithosphere__increment_of_overlying_pressure", at="node")
    return Flexure1D(grid)


FLEXURE_KEYWORDS = ("eet", "youngs", "rho_mantle", "rho_water", "gravity")


def pytest_generate_tests(metafunc):
    if "flexure_keyword" in metafunc.fixturenames:
        metafunc.parametrize("flexure_keyword", FLEXURE_KEYWORDS)

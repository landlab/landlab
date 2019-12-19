import pytest

from landlab import RasterModelGrid
from landlab.components.plant_competition_ca.plant_competition_ca import VegCA


@pytest.fixture
def ca_veg():
    grid = RasterModelGrid((20, 20), xy_spacing=10e0)
    return VegCA(grid)

import pytest

from landlab import RasterModelGrid
from landlab.components.plant_competition_ca.plant_competition_ca import VegCA


@pytest.fixture
def ca_veg():
    grid = RasterModelGrid((20, 20), xy_spacing=10e0)
    grid.add_zeros("vegetation__cumulative_water_stress", at="cell")
    grid.add_zeros("vegetation__plant_functional_type", at="cell", dtype=int)
    return VegCA(grid)

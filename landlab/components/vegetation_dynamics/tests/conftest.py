import pytest

from landlab import RasterModelGrid
from landlab.components.vegetation_dynamics.vegetation_dynamics import Vegetation


@pytest.fixture
def veg():
    grid = RasterModelGrid((20, 20), xy_spacing=10e0)
    grid.add_zeros("vegetation__plant_functional_type", at="cell", dtype=int)
    return Vegetation(grid)

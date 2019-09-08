import pytest

from landlab import RasterModelGrid
from landlab.components.vegetation_dynamics.vegetation_dynamics import Vegetation


@pytest.fixture
def veg():
    grid = RasterModelGrid((20, 20), xy_spacing=10e0)
    grid.add_zeros("vegetation__plant_functional_type", at="cell", dtype=int)
    grid.add_zeros("surface__evapotranspiration", at="cell", dtype=float)
    grid.add_zeros("vegetation__water_stress", at="cell", dtype=float)
    grid.add_zeros("surface__potential_evapotranspiration_rate", at="cell", dtype=float)
    grid.add_zeros(
        "surface__potential_evapotranspiration_30day_mean", at="cell", dtype=float
    )

    return Vegetation(grid)

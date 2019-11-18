import pytest

from landlab import RasterModelGrid
from landlab.components.pet.potential_evapotranspiration_field import (
    PotentialEvapotranspiration,
)


@pytest.fixture
def pet():
    grid = RasterModelGrid((20, 20), xy_spacing=10e0)
    return PotentialEvapotranspiration(grid)

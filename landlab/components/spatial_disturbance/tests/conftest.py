
import pytest

from landlab import RasterModelGrid
from landlab.components.spatial_disturbance.spatial_disturbance import SpatialDisturbance


@pytest.fixture
def sd():
    grid = RasterModelGrid((20, 20), xy_spacing=10e0)
    grid.add_zeros("vegetation__plant_functional_type", at="cell", dtype=int)
    return SpatialDisturbance(grid)

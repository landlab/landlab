import pytest

from landlab import RasterModelGrid
from landlab.components.resource_redistribution.resource_redistribution import ResourceRedistribution


@pytest.fixture
def rr():
    grid = RasterModelGrid((20, 20), xy_spacing=10e0)
    grid.add_zeros("vegetation__plant_functional_type", at="cell", dtype=int)
    grid.add_zeros("soil__resources", at="cell", dtype=float)
    return ResourceRedistribution(grid)

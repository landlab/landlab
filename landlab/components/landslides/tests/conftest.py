import pytest

from landlab import RasterModelGrid
from landlab.components import LandslideProbability


@pytest.fixture
def ls_prob():
    grid = RasterModelGrid((20, 20), xy_spacing=10e0)
    grid.add_zeros("topographic__slope", at="node", dtype=float)
    return LandslideProbability(grid)

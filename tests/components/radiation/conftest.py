import pytest

from landlab import RasterModelGrid
from landlab.components.radiation.radiation import Radiation


@pytest.fixture
def rad():
    grid = RasterModelGrid((20, 20), xy_spacing=10e0)
    grid.add_zeros("topographic__elevation", at="node")
    return Radiation(grid)

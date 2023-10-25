import pytest

from landlab import RasterModelGrid
from landlab.components import RiverBedDynamics


@pytest.fixture
def r_b_d():
    grid = RasterModelGrid((32, 240), xy_spacing=25)
    grid.add_zeros("bed_surface__grain_size_distribution_location", at="node")
    grid.add_zeros("surface_water__depth", at="link")
    grid.add_zeros("surface_water__velocity", at="link")
    grid.add_zeros("surface_water__velocity_previous_time", at="link")
    grid.add_zeros("topographic__elevation", at="node")
    return RiverBedDynamics(grid)

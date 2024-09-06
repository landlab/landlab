import pytest

from landlab import RasterModelGrid
from landlab.components import RiverBedDynamics


@pytest.fixture
def rbd():
    grid = RasterModelGrid((5, 5), xy_spacing=100)

    grid.at_node["topographic__elevation"] = [
        [1.07, 1.08, 1.09, 1.09, 1.09],
        [1.06, 1.07, 1.08, 1.09, 1.09],
        [1.00, 1.03, 1.07, 1.08, 1.09],
        [1.06, 1.07, 1.08, 1.09, 1.09],
        [1.07, 1.08, 1.09, 1.09, 1.09],
    ]

    grid.set_watershed_boundary_condition(grid.at_node["topographic__elevation"])

    grid.add_zeros("bed_surface__grain_size_distribution_location", at="node")
    grid.add_zeros("surface_water__depth", at="link")
    grid.add_zeros("surface_water__velocity", at="link")
    grid["link"]["surface_water__velocity"][20] = 10
    grid.add_zeros("surface_water__velocity_previous_time", at="link")

    gsd = [[128, 100], [64, 90], [32, 80], [16, 50], [8, 20], [4, 10], [2, 1], [1, 0]]

    return RiverBedDynamics(grid, gsd=gsd)

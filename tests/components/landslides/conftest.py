import pytest

from landlab import RasterModelGrid
from landlab.components import LandslideProbability


@pytest.fixture
def ls_prob():
    grid = RasterModelGrid((20, 20), xy_spacing=10e0)
    grid.add_zeros("topographic__slope", at="node", dtype=float)
    grid.add_zeros("topographic__specific_contributing_area", at="node")
    grid.add_zeros("soil__transmissivity", at="node")
    grid.add_zeros("soil__saturated_hydraulic_conductivity", at="node")
    grid.add_zeros("soil__mode_total_cohesion", at="node")
    grid.add_zeros("soil__minimum_total_cohesion", at="node")
    grid.add_zeros("soil__maximum_total_cohesion", at="node")
    grid.add_zeros("soil__internal_friction_angle", at="node")
    grid.add_zeros("soil__density", at="node")
    grid.add_zeros("soil__thickness", at="node")

    return LandslideProbability(grid)

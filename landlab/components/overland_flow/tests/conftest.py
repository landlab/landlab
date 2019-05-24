import pytest

from landlab import RasterModelGrid
from landlab.components.overland_flow import (
    KinwaveOverlandFlowModel,
    OverlandFlow,
    OverlandFlowBates,
)


@pytest.fixture
def deAlm():
    grid = RasterModelGrid((32, 240), xy_spacing=25)
    grid.add_zeros("surface_water__depth", at="node")
    grid.add_zeros("topographic__elevation", at="node")
    grid.add_zeros("surface_water__discharge", at="link")
    return OverlandFlow(grid, mannings_n=0.01, h_init=0.001)


@pytest.fixture
def kin_wave_of():
    grid = RasterModelGrid((10, 10), xy_spacing=0.5)
    grid.add_zeros("topographic__elevation", at="node", dtype=float)
    grid.add_zeros("topographic__gradient", at="node")

    return KinwaveOverlandFlowModel(grid)


@pytest.fixture
def bates():
    grid = RasterModelGrid((32, 240), xy_spacing=25)
    grid.add_zeros("surface_water__depth", at="node")
    grid.add_zeros("topographic__elevation", at="node")
    return OverlandFlowBates(grid, mannings_n=0.01, h_init=0.001)

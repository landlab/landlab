import pytest

from landlab import RasterModelGrid
from landlab.components.gflex.flexure import gFlex

pytest.importorskip("gflex")


@pytest.fixture
def grid():
    mg = RasterModelGrid((20, 20), xy_spacing=25000.0)
    mg.add_zeros("surface_load__stress", at="node")
    return mg


@pytest.fixture
def grid_with_topo():
    mg = RasterModelGrid((20, 20), xy_spacing=25000.0)
    mg.add_zeros("surface_load__stress", at="node")
    mg.add_zeros("topographic__elevation", at="node")
    return mg


@pytest.fixture
def gf(grid):
    return gFlex(grid, quiet=True)

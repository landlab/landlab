import pytest

from landlab import RasterModelGrid
from landlab.components.gflex.flexure import gFlex

pytest.importorskip("gflex")


@pytest.fixture
def grid():
    mg = RasterModelGrid((20, 20), xy_spacing=25000.0)
    mg.add_zeros("load__normal_component_of_stress", at="node")
    return mg


@pytest.fixture
def gf(grid):
    return gFlex(grid, quiet=True)

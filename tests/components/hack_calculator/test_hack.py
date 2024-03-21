import pytest

from landlab import RasterModelGrid
from landlab.components import FlowAccumulator
from landlab.components import HackCalculator


@pytest.fixture()
def simple_hack_test_grid():
    mg = RasterModelGrid((10, 10))
    z = mg.add_zeros("topographic__elevation", at="node")
    z += mg.x_of_node + mg.y_of_node
    fa = FlowAccumulator(mg)
    fa.run_one_step()
    return mg


def test_no_full_hack_dataframe_specified(simple_hack_test_grid):
    mg = simple_hack_test_grid
    hc = HackCalculator(mg)
    with pytest.raises(NotImplementedError):
        hc.full_hack_dataframe


def test_no_full_hack_dataframe_built(simple_hack_test_grid):
    mg = simple_hack_test_grid
    hc = HackCalculator(mg, save_full_df=True)
    with pytest.raises(RuntimeError):
        hc.full_hack_dataframe


def test_no_coefficient_df_built(simple_hack_test_grid):
    mg = simple_hack_test_grid
    hc = HackCalculator(mg)
    with pytest.raises(RuntimeError):
        hc.hack_coefficient_dataframe

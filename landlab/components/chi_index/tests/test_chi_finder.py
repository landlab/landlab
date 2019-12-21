import pytest

from landlab import HexModelGrid, RasterModelGrid
from landlab.components import ChiFinder, FlowAccumulator


def test_route_to_multiple_error_raised():
    mg = RasterModelGrid((10, 10))
    z = mg.add_zeros("node", "topographic__elevation")
    z += mg.x_of_node + mg.y_of_node
    fa = FlowAccumulator(mg, flow_director="MFD")
    fa.run_one_step()

    with pytest.raises(NotImplementedError):
        ChiFinder(mg, min_drainage_area=1.0, reference_concavity=1.0)


def test_bad_reference_area():
    mg = RasterModelGrid((10, 10))
    z = mg.add_zeros("node", "topographic__elevation")
    z += mg.x_of_node + mg.y_of_node
    fa = FlowAccumulator(mg)
    fa.run_one_step()

    with pytest.raises(ValueError):
        ChiFinder(mg, min_drainage_area=1.0, reference_area=-1)


def test_functions_with_Hex():
    mg = HexModelGrid(10, 10)
    z = mg.add_zeros("node", "topographic__elevation")
    z += mg.x_of_node + mg.y_of_node
    fa = FlowAccumulator(mg)
    fa.run_one_step()

    ch = ChiFinder(mg, min_drainage_area=1.0, reference_concavity=1.0)
    ch.calculate_chi()

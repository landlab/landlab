import matplotlib.pyplot as plt

from landlab import RasterModelGrid
from landlab.components import FlowAccumulator
from landlab.plot.drainage_plot import drainage_plot


def make_grid():
    mg = RasterModelGrid((10, 10), xy_spacing=(1, 1))
    mg.add_field(
        "topographic__elevation", mg.node_x ** 2 + mg.node_y ** 2 + mg.node_y, at="node"
    )
    return mg


def test_steepest():
    mg = make_grid()
    fa = FlowAccumulator(mg)
    fa.run_one_step()
    plt.figure()
    drainage_plot(mg)


def test_mfd():
    mg = make_grid()
    fa = FlowAccumulator(mg, flow_director="MFD")
    fa.run_one_step()
    plt.figure()
    drainage_plot(mg)

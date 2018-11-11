import numpy as np
from matplotlib.pyplot import show
from six.moves import range

from landlab import VoronoiDelaunayGrid
from landlab.components import FlowAccumulator, StreamPowerEroder
from landlab.plot.imshow import imshow_node_grid

nnodes = 10000

x, y = np.random.rand(nnodes), np.random.rand(nnodes)
mg = VoronoiDelaunayGrid(x, y)

z = mg.add_field(
    "node", "topographic__elevation", np.random.rand(nnodes) / 10000., copy=False
)

fr = FlowAccumulator(mg)
spe = StreamPowerEroder(mg, "drive_sp_params_voronoi.txt")

for i in range(100):
    z[mg.core_nodes] += 0.01
    fr.run_one_step()
    spe.erode(mg, 1.)

imshow_node_grid(mg, "topographic__elevation")

show()

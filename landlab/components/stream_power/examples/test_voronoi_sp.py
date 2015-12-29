from six.moves import range

from landlab import VoronoiDelaunayGrid
from landlab.components.flow_routing import FlowRouter
from landlab.components.stream_power import StreamPowerEroder
from landlab.plot.imshow import imshow_node_grid
import numpy as np
from matplotlib.pyplot import figure, show

nnodes = 10000

x, y = np.random.rand(nnodes), np.random.rand(nnodes)
mg = VoronoiDelaunayGrid(x,y)

z = mg.add_field('node', 'topographic__elevation', np.random.rand(nnodes)/10000., copy=False)

fr = FlowRouter(mg)
spe = StreamPowerEroder(mg, 'drive_sp_params_voronoi.txt')

for i in range(100):
    z[mg.core_nodes] += 0.01
    fr.route_flow()
    spe.erode(mg, 1.)

imshow_node_grid(mg, 'topographic__elevation')

show()

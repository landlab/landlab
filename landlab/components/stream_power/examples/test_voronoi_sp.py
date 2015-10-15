from landlab import VoronoiDelaunayGrid  # , RasterModelGrid
from landlab.components.flow_routing.route_flow_dn import FlowRouter
from landlab.components.stream_power.stream_power import StreamPowerEroder
import numpy as np

x, y = np.random.rand(50), np.random.rand(50)
mg = VoronoiDelaunayGrid(x,y)
#mg = RasterModelGrid(4,5)

mg.add_field('node', 'topographic__elevation', mg.node_x, copy=True)

fr = FlowRouter(mg)
spe = StreamPowerEroder(mg, 'drive_sp_params.txt')

for i in xrange(100):
    fr.route_flow()
    spe.erode(1.)

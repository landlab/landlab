from six.moves import range
import numpy as np
from matplotlib.pyplot import show, plot, figure
from landlab import RasterModelGrid, CLOSED_BOUNDARY, imshow_grid_at_node
from landlab.components import FlowRouter, SedDepEroder

mg = RasterModelGrid((30, 3), 200.)  # ((10, 3), 200.)
for edge in (mg.nodes_at_left_edge, mg.nodes_at_top_edge,
             mg.nodes_at_right_edge):
    mg.status_at_node[edge] = CLOSED_BOUNDARY

z = mg.add_zeros('node', 'topographic__elevation')
th = mg.add_zeros('node', 'channel_sediment__depth')
th += 0.0007

fr = FlowRouter(mg)
sde = SedDepEroder(mg, K_sp=1.e-5,
                   sed_dependency_type='almost_parabolic',
                   Qc='power_law', K_t=1.e-5)

z[:] = mg.node_y/10000.

initz = z.copy()

dt = 100.
up = 0.0005

for i in range(12000):
    z[mg.core_nodes] += up*dt
    fr.run_one_step()
    sde.run_one_step(dt)
    if i%1000 == 0:
        plot(z.reshape((30, 3))[:-1, 1])

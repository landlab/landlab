from six.moves import range
import numpy as np
from matplotlib.pyplot import show, plot, figure
from landlab import RasterModelGrid, CLOSED_BOUNDARY, imshow_grid_at_node
from landlab.components import FlowRouter, SedDepEroder, DepressionFinderAndRouter

mg = RasterModelGrid((10, 3), 200.)
for edge in (mg.nodes_at_left_edge, mg.nodes_at_top_edge,
             mg.nodes_at_right_edge):
    mg.status_at_node[edge] = CLOSED_BOUNDARY

z = mg.add_zeros('node', 'topographic__elevation')
th = mg.add_zeros('node', 'channel_sediment__depth')
th += 0.001

fr = FlowRouter(mg)
pit = DepressionFinderAndRouter(mg)
sde = SedDepEroder(mg, K_sp=1.e-5,
                   sed_dependency_type='almost_parabolic',
                   Qc='power_law', K_t=1.e-5)

z[:] = mg.node_y/10000.
# make a rock dam:
#z[16] -= 0.08 #0.065
#z[13] -= 0.05
#z[10] -= 0.022

dt = 100.

#myflood = np.arange(1, 10, 3)  # 1st 3 nodes, incl boundary, are flooded
# or as bool
myflood = mg.zeros('node', dtype=bool)
myflood[mg.node_y < 500.] = True

initz = z.copy()

for i in range(300):
    fr.run_one_step()
    pit.map_depressions()
    sde.run_one_step(dt, flooded_nodes=myflood)
    if i%10 == 0:
        figure(1)
        plot(z.reshape((10, 3))[:-1, 1])
        figure(2)
        plot(th.reshape((10, 3))[:-1, 1])
        figure(3)
        plot(pit.depression_depth.reshape((10, 3))[:-1, 1])

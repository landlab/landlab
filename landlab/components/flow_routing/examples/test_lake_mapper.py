from landlab import RasterModelGrid
from landlab.plot.imshow import imshow_node_grid
from landlab.components.flow_routing import (FlowRouter,
                                             DepressionFinderAndRouter)
from matplotlib.pyplot import figure, plot, show
import numpy as np

nx, ny = 50, 50
mg = RasterModelGrid(nx,ny,1.)

# make a "crater"
x_distance_from_center = mg.node_x-mg.node_x.mean()
y_distance_from_center = mg.node_y-mg.node_y.mean()
x_distance_from_edge = np.amin(np.vstack((mg.node_x.max()-mg.node_x,
                                          mg.node_x-mg.node_x.min())), axis=0)
y_distance_from_edge = np.amin(np.vstack((mg.node_y.max()-mg.node_y,
                                          mg.node_y-mg.node_y.min())), axis=0)
# make the "hole"
hole_elev = np.sqrt(x_distance_from_center**2 + y_distance_from_center**2)
# make the rim:
rim_elev = np.sqrt(x_distance_from_edge**2 + y_distance_from_edge**2)
# assemble
z = np.amin(np.vstack((hole_elev, rim_elev)), axis=0)
z += np.random.rand(nx*ny)/1000.

mg.add_field('node', 'topographic__elevation', z, copy=False)

fr = FlowRouter(mg)
lf = DepressionFinderAndRouter(mg)

fr.route_flow()

figure('old drainage area')
imshow_node_grid(mg, 'drainage_area')

lf.map_depressions(pits=mg.at_node['flow__sink_flag'])

figure('depression depth')
imshow_node_grid(mg, 'depression__depth')

figure('new drainage area')
imshow_node_grid(mg, 'drainage_area')

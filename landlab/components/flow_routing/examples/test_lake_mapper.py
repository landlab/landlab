from landlab import RasterModelGrid
from landlab.plot.imshow import imshow_node_grid
from landlab.components.flow_routing.route_flow_dn import FlowRouter
from landlab.components.flow_routing.lake_mapper import DepressionFinderAndRouter
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

# This configuration can become a unit test:

# z = np.array([  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
#                 0.0,  2.0,  2.0,  2.0,  2.0,  2.0,  0.0,
#                 0.0,  2.0,  1.6,  1.5,  1.6,  2.0,  0.0,
#                 0.0,  2.0,  1.7,  1.6,  1.7,  2.0,  0.0,
#                 0.0,  2.0,  1.8,  2.0,  2.0,  2.0,  0.0,
#                 0.0,  1.0,  0.6,  1.0,  1.0,  1.0,  0.0,
#                 0.0,  0.0, -0.5,  0.0,  0.0,  0.0,  0.0])

r_old = np.array([  0,  1,  2,  3,  4,  5,  6,
                    7,  1,  2,  3,  4,  5, 13,
                    14, 14, 17, 17, 17, 20, 20,
                    21, 21, 17, 17, 17, 27, 27,
                    28, 28, 37, 37, 39, 34, 34,
                    35, 44, 44, 44, 46, 47, 41,
                    42, 43, 44, 45, 46, 47, 48])

r_new = np.array([  0,  1,  2,  3,  4,  5,  6,
                    7,  1,  2,  3,  4,  5, 13,
                    14, 14, 17, 17, 17, 20, 20,
                    21, 21, 30, 30, 17, 27, 27,
                    28, 28, 37, 37, 39, 34, 34,
                    35, 44, 44, 44, 46, 47, 41,
                    42, 43, 44, 45, 46, 47, 48])

mg.add_field('node', 'topographic__elevation', z, copy=False)

fr = FlowRouter(mg)
lf = DepressionFinderAndRouter(mg)

fr.route_flow()

figure('old drainage area')
imshow_node_grid(mg, 'drainage_area')

lf.map_depressions()

figure('depression depth')
imshow_node_grid(mg, 'depression__depth')

figure('new drainage area')
imshow_node_grid(mg, 'drainage_area')

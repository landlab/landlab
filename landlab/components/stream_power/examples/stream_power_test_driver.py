"""
This is not yet working. Hopefully updating the D8 flow routing methods in the grid will resolve this.
"""

import numpy as np
from pylab import show, imshow, colorbar, plot
from landlab.grid.raster import RasterModelGrid
from landlab.components.flow_routing import flow_routing_D8
from landlab.components.flow_accum import flow_accumulation
from landlab.components.stream_power import stream_power
from landlab.plot.imshow import imshow_grid
reload(stream_power)

class data(object):
    '''
        This is where all the whole-grid data lives, as arrays over the various elements of the grid.
        '''
    #Data goes here!!!
    def __init__(self, grid):
        self.elev = grid.create_node_dvector() #some data

iterations = 10
tstep = 1.

#Make grid, set the elevs
mg = RasterModelGrid()
mg.initialize(5,10,5.)
mg.set_inactive_boundaries(True, True, True, True)
vectors = data(mg)

#A surface dipping right, plus random noise:
loading_vector = np.linspace(10,1,num=10)
vectors.elev = np.concatenate((loading_vector,loading_vector,loading_vector,loading_vector,loading_vector))
vectors.elev += np.random.random_sample(vectors.elev.shape)/10.

print vectors.elev.shape
print type(vectors.elev)
print vectors.elev.reshape((5,10))
print mg.node_status.reshape((5,10))

#Initialize the modules:
dirs = flow_routing_D8.RouteFlowD8(len(vectors.elev))
accum = flow_accumulation.AccumFlow(mg, vectors)
sp = stream_power.StreamPower(mg, vectors, tstep)

#Run the modules:
for i in xrange(iterations):
    vectors.flowdirs = dirs.calc_flowdirs(mg, vectors.elev)
    vectors.flowacc = accum.calc_flowacc(mg, vectors)
    vectors.elev = sp.stream_power_erosion(mg, vectors)
    
print vectors.elev.reshape((5,10))
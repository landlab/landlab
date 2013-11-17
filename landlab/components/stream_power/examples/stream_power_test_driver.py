"""
This is not yet working. Hopefully updating the D8 flow routing methods in the grid will resolve this.
"""

import numpy as np
from pylab import show, imshow, colorbar, plot
import landlab.grid.raster as raster
from landlab.components.flow_routing import flow_routing_D8
from landlab.components.flow_accum import flow_accumulation
from landlab.components.stream_power import stream_power
from landlab.plot.imshow import imshow_grid
from pylab import show
from copy import copy
reload(flow_routing_D8)
reload(flow_accumulation)
reload(stream_power)
reload(raster)

class data(object):
    '''
        This is where all the whole-grid data lives, as arrays over the various elements of the grid.
        '''
    #Data goes here!!!
    def __init__(self, grid):
        self.elev = grid.zeros() #some data

iterations = 10
tstep = 1.

#Make grid, set the elevs
mg = raster.RasterModelGrid(50, 50, 1.)
mg.set_inactive_boundaries(True, True, True, True)
vectors = data(mg)

#A surface dipping right, plus random noise:
loading_vector = np.linspace(50,1,num=50)
vectors.elev = np.repeat(loading_vector, 50)
vectors.elev += np.random.random_sample(vectors.elev.shape)/10.
vectors.init_elev = copy(vectors.elev)

print vectors.elev.shape
#print vectors.elev.reshape((5,10))

print vectors.elev.shape
print type(vectors.elev)
#print vectors.elev.reshape((5,10))
#print mg.node_status.reshape((5,10))

#Initialize the modules:
dirs = flow_routing_D8.RouteFlowD8(len(vectors.elev))
accum = flow_accumulation.AccumFlow(mg, vectors)
sp = stream_power.StreamPower(mg, vectors, tstep)

#print vectors.elev.reshape((5,10))

#Run the modules:
for i in xrange(iterations):
    print i
    vectors.flowdirs, vectors.node_max_gradients = dirs.calc_flowdirs(mg, vectors.elev)
    try:
        print type(vectors.elev)
    except:
        print 'broke after flowdirs'
    vectors.flowacc = accum.calc_flowacc(mg, vectors)
    try:
        print type(vectors.elev)
    except:
        print 'broke after flowacc'
    vectors.elev = sp.stream_power_erosion(mg, vectors)
    try:
        print type(vectors.elev)
        print vectors.elev.shape
    except:
        print 'broke after sp'
    
#print vectors.elev.reshape((5,10))

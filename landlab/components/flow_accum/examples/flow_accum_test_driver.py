import numpy as np
from pylab import show, imshow, colorbar, plot
import landlab.grid.raster as raster
from landlab.components.flow_routing import flow_routing_D8
from landlab.components.flow_accum import flow_accumulation
from landlab.plot.imshow import imshow_grid
reload(flow_routing_D8)
reload(flow_accumulation)
reload(raster)

class data(object):
    '''
        This is where all the whole-grid data lives, as arrays over the various elements of the grid.
        '''
    #Data goes here!!!
    def __init__(self, grid):
        self.elev = grid.zeros(centering='node') #some data

#Make grid, set the elevs
mg = raster.RasterModelGrid(5, 5, 1.)
mg.set_inactive_boundaries(True, True, False, True)
vectors = data(mg)

vectors.elev = np.array([10.,10.,0.,10.,10.,10.,8.,8.,9.,10.,10.,4.,6.,7.,10.,10.,1.,3.,8.,10.,10.,0.,10.,10.,10.])

print vectors.elev.reshape((5,5))
print np.array(range(25)).reshape((5,5))
#print mg.node_status.reshape((5,5))

#imshow_grid(mg,vectors.elev)
#plot()

#Route the flow
dirs = flow_routing_D8.RouteFlowD8(len(vectors.elev))
vectors.flowdirs, vectors.node_max_gradients = dirs.calc_flowdirs(mg, vectors.elev)

#print vectors.flowdirs.shape
#print vectors.flowdirs.reshape((5,5))

#Accumulate the flow
accum = flow_accumulation.AccumFlow(mg, vectors)
#vectors.flowacc = accum.calc_flowacc(mg, vectors)
accum.calc_flowacc(mg, vectors)

print vectors.flowacc.shape
print vectors.flowacc.reshape((5,5))

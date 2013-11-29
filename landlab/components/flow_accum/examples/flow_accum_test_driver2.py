import numpy as np
from pylab import show, imshow, colorbar, plot
from landlab import RasterModelGrid
from landlab.components.flow_routing.flow_routing_D8 import RouteFlowD8
from landlab.components.flow_accum.flow_accumulation2 import AccumFlow
from landlab.plot.imshow import imshow_grid
from landlab.components.dem_support.dem_boundary_conditions import WatershedBoundaryConditions
from random import uniform
#reload(flow_routing_D8)
#reload(flow_accumulation)
#reload(raster)

def main():
    nr = 5
    nc = 6
    nnodes = nr*nc
    dx=3
    #instantiate grid
    rg = RasterModelGrid(nr, nc, dx)
    #rg.set_inactive_boundaries(False, False, True, True)
    
    nodata_val=-1
    z  = nodata_val*np.ones( nnodes )    
    #set-up interior elevations with random numbers
    #for i in range(0, nnodes):
    #    if rg.is_interior(i):
    #        elevations[i]=random.random_sample()
    
    #set-up with prescribed elevations to test drainage area calcualtion
    helper = [7,8,9,10,13,14,15,16]
    for i in xrange(0, len(helper)):
        #print 'helper[i]', helper[i]
        z[helper[i]]=2+uniform(-0.5,0.5)
    helper = [19,20,21,22]
    for i in xrange(0, len(helper)):
        z[helper[i]]=3+uniform(-0.5,0.5)
        
    z[7]=1
    
    bc=WatershedBoundaryConditions()
    bc.set_bc_find_outlet(rg, z, nodata_val)
    
    #instantiate variable of type RouteFlowD8 Class
    flow_router = RouteFlowD8(len(z))
    #initial flow direction
    flowdirs, max_slopes = flow_router.calc_flowdirs(rg,z)
    #insantiate variable of type AccumFlow Class
    accumulator = AccumFlow(rg)
    #initial flow accumulation
    drain_area = accumulator.calc_flowacc(rg, z, flowdirs)
    
    print "elevations ", rg.node_vector_to_raster(z)
    print "flowdirs ", rg.node_vector_to_raster(flowdirs)
    print "drain_area ", rg.node_vector_to_raster(drain_area)
    
    
if __name__ == '__main__':
    main()

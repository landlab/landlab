import numpy as np
from pylab import *
from landlab import RasterModelGrid, CLOSED_BOUNDARY
from landlab.plot.imshow import imshow_grid
from landlab.components.dem_support import WatershedBoundaryConditions
from random import uniform
from landlab.components.simple_power_law_incision import PowerLawIncision
from landlab.components.flow_routing import RouteFlowD8
from landlab.components.flow_accum import AccumFlow
import matplotlib.pyplot as plt


def main():
    nr = 5
    nc = 6
    nnodes = nr*nc
    dx=1
    #instantiate grid
    rg = RasterModelGrid(nr, nc, dx)
    #rg.set_inactive_boundaries(False, False, True, True)

    nodata_val=-1
    z  = nodata_val*np.ones( nnodes )
    #set-up interior elevations with random numbers
    #for i in range(0, nnodes):
    #    if rg.is_interior(i):
    #        elevations[i]=random.random_sample()

    #set-up elevations
    helper = [7,8,9,10,13,14,15,16]
    for i in xrange(0, len(helper)):
        #print 'helper[i]', helper[i]
        z[helper[i]]=2+uniform(-0.5,0.5)
    helper = [19,20,21,22]
    for i in xrange(0, len(helper)):
        z[helper[i]]=3+uniform(-0.5,0.5)

    z[7]=1

    #set up boundary conditions
    bc=WatershedBoundaryConditions()
    outlet_loc = bc.set_bc_find_outlet(rg, z, nodata_val)
    zoutlet=z[outlet_loc]

    #some helper and parameter values
    total_run_time = 500000 #yr
    uplift_rate = 0.001 #m/yr
    rain_rate = 1 #m/yr
    storm_duration = 50 #years

    elapsed_time = 0 #years

    #set up fluvial incision component
    incisor = PowerLawIncision('input_file.txt',rg)

    #print "elevations before ", rg.node_vector_to_raster(z)
    #interior_nodes are the nodes on which you will be operating
    interior_nodes = np.where(rg.status_at_node != CLOSED_BOUNDARY)[0]

    while elapsed_time < total_run_time:
        #uplift the landscape
        z[interior_nodes] = z[interior_nodes]+uplift_rate * storm_duration
        z[outlet_loc]=zoutlet
        #erode the landscape
        z = incisor.run_one_storm(rg,z,rain_rate,storm_duration)
        #update the time
        elapsed_time = elapsed_time+storm_duration

    #below purely for plotting reasons
    #instantiate variable of type RouteFlowD8 Class
    flow_router = RouteFlowD8(len(z))
    #initial flow direction
    flowdirs, max_slopes = flow_router.calc_flowdirs(rg,z)
    #insantiate variable of type AccumFlow Class
    accumulator = AccumFlow(rg)
    #initial flow accumulation
    drain_area = accumulator.calc_flowacc(z, flowdirs)

    #m,b = polyfit(log10(drain_area[interior_nodes]), log10(max_slopes[interior_nodes]), 1)
    z[interior_nodes] = z[interior_nodes]+uplift_rate * storm_duration
    z[outlet_loc]=zoutlet

    plt.loglog(np.array(drain_area),np.array(max_slopes),'ro',)
    plt.show()
    imshow_grid(rg,z, values_at='node')
    #print "elevations after", rg.node_vector_to_raster(z)
    #print "drainage area after", rg.node_vector_to_raster(drain_area)
    #print "slope after", rg.node_vector_to_raster(max_slopes)
    #print "flowdirs after", rg.node_vector_to_raster(flowdirs)


if __name__ == '__main__':
    main()

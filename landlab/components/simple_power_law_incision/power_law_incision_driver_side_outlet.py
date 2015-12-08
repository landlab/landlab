from __future__ import print_function

import numpy as np
import pylab
from landlab import RasterModelGrid, CLOSED_BOUNDARY
from random import uniform
from landlab.components.simple_power_law_incision import PowerLawIncision
from landlab.components.flow_routing import RouteFlowD8
from landlab.components.flow_accum import AccumFlow
import matplotlib.pyplot as plt


def main():
    nr = 50
    nc = 60
    nnodes = nr*nc
    dx=10
    #instantiate grid
    rg = RasterModelGrid(nr, nc, dx)
    rg.set_inactive_boundaries(False, True, True, True)

    z  = np.zeros( nnodes )
    #set-up interior elevations with random numbers
    for i in range(0, nnodes):
        if rg.is_interior(i):
            z[i]=2+uniform(-0.5,0.5)

    #some helper and parameter values
    total_run_time = 10000 #yr
    one_twentieth_time = total_run_time/20
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
        #erode the landscape
        z = incisor.run_one_storm(rg,z,rain_rate,storm_duration)
        #uplift the landscape
        z[interior_nodes] = z[interior_nodes]+uplift_rate * storm_duration
        #update the time
        elapsed_time = elapsed_time+storm_duration
        #print "elapsed_time", elapsed_time
        if elapsed_time%one_twentieth_time == 0:
            print("elapsed time",elapsed_time)
            elev_raster = rg.node_vector_to_raster(z,True)
            # Plot topography
            pylab.figure(22)
            im = pylab.imshow(elev_raster, cmap=pylab.cm.RdBu, extent=[0, nc*rg.dx, 0, nr*rg.dx])
            cb = pylab.colorbar(im)
            cb.set_label('Elevation (m)', fontsize=12)
            pylab.title('Topography')
            pylab.show()

    #below purely for plotting reasons
    #instantiate variable of type RouteFlowD8 Class
    flow_router = RouteFlowD8(len(z))
    #initial flow direction
    flowdirs, max_slopes = flow_router.calc_flowdirs(rg,z)
    #insantiate variable of type AccumFlow Class
    accumulator = AccumFlow(rg)
    #initial flow accumulation
    drain_area = accumulator.calc_flowacc(z, flowdirs)

    plt.loglog(np.array(drain_area),np.array(max_slopes),'ro',)
    plt.xlabel('drainage area, m')
    plt.ylabel('surface slope')
    plt.show()

    elev_raster = rg.node_vector_to_raster(z,True)
    # Plot topography
    pylab.figure(22)
    im = pylab.imshow(elev_raster, cmap=pylab.cm.RdBu, extent=[0, nc*rg.dx, 0, nr*rg.dx])
    cb = pylab.colorbar(im)
    cb.set_label('Elevation (m)', fontsize=12)
    pylab.title('Topography')

    pylab.show()


if __name__ == '__main__':
    main()

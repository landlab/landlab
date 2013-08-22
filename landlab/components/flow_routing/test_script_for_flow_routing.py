#! /usr/env/python
"""
just a little script for testing the d8 flow routing class 
and drainage area calculation 
"""

from landlab import RasterModelGrid
import flow_routing_D8
from landlab.examples.drainageArea import CalcDrainageArea
from numpy import *
from pylab import plot, draw, show, contour, imshow, colorbar

def main():
    nr = 5
    nc = 6
    ncells = nr*nc
    dx=1
    #instantiate grid
    rg = RasterModelGrid(nr, nc, dx)
    rg.set_inactive_boundaries(False, False, True, True)
    

    elevations  = zeros( ncells )    
    #set-up interior elevations with random numbers
    #for i in range(0, ncells):
    #    if rg.is_interior(i):
    #        elevations[i]=random.random_sample()
    
    #set-up with prescribed elevations to test drainage area calcualtion
    helper = [7,8,9,10,13,14,15,16]
    for i in range(0, len(helper)):
        #print 'helper[i]', helper[i]
        elevations[helper[i]]=2
    helper = [19,20,21,22]
    for i in range(0, len(helper)):
        elevations[helper[i]]=3
        
    elevations[7]=1
    
    #tried making a pit, drainage area algorithm doesn't crash, but it doesn't
    #work perfectly either.    
    #elevations[15]=-50
    
    #printing elevations for debugging purposes
    #print 'elevation vector' 
    #print elevations

    #instantiate flow routing variable
    flow = flow_routing_D8.RouteFlowD8(ncells)
    #calculate flow directions
    flow_directions = flow.calc_flowdirs(rg, elevations) 
    fd_raster = rg.node_vector_to_raster(flow_directions,True)   
    #printing flow directions for debugging purposes
    #print 'flow direction vector'
    #print flow_directions

    
    #instantiate drainage area variable
    da_calculator = CalcDrainageArea(ncells)
    #calculate drainage area
    drain_area = da_calculator.calc_DA(rg, flow_directions)
    #printing drainage area for debugging purposes
    #print 'drainage area vector'
    #print drain_area
    
    #printing, from Dan's craters code
    elev_raster = rg.node_vector_to_raster(elevations,True)
    #contour(elev_raster)
    #flipped_elev_raster = numpy.empty_like(elev_raster)
    #for i in range(0,nr):
    #    flipped_elev_raster[i,:] = elev_raster[(nr-i-1),:]
    #imshow(flipped_elev_raster)
    print 'elevation raster'
    print elev_raster
    #imshow(elev_raster)
    #colorbar()
    #show()
    
    print 'flow direction raster'
    print fd_raster
    
    da_raster = rg.node_vector_to_raster(drain_area)
    #contour(elev_raster)
    #flipped_elev_raster = numpy.empty_like(elev_raster)
    #for i in range(0,nr):
    #    flipped_elev_raster[i,:] = elev_raster[(nr-i-1),:]
    #imshow(flipped_elev_raster)
    #imshow(da_raster)
    #colorbar()
    #show()
    
    
    ##test of other code bit in raster model grid
    #[ms, ma] = rg.calculate_max_gradient_across_node(elevations,12)
    #print 'max slope', ms
    #print 'max angle', ma


if __name__ == '__main__':
    main()

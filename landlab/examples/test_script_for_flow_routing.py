#! /usr/env/python
"""
just a little script for testing the d8 flow routing class 
"""

from landlab.model_grid import RasterModelGrid
from landlab.examples.flowRoutingD8 import RouteFlowD8
from numpy import *

def main():
    nr = 5
    nc = 5
    ncells = nr*nc
    dx=1
    #instantiate grid
    rg = RasterModelGrid(nr, nc, dx)
    
    #set-up interior elevations with random numbers
    elevations  = zeros( ncells )
    for i in range(0, ncells):
        if rg.is_interior(i):
            elevations[i]=random.random_sample()
    
    #printing elevations for debugging purposes
    print 'elevation vector' 
    print elevations

    #instantiate flow routing variable
    flow = RouteFlowD8(ncells)
    #calculate flow directions
    flow_directions = flow.calc_flowdirs(rg, elevations)    
    #printing flow directions for debugging purposes
    print 'flow direction vector'
    print flow_directions
    
    #test of other code bit in raster model grid
    [ms, ma] = rg.calculate_max_gradient_across_node(elevations,12)
    print 'max slope', ms
    print 'max angle', ma


if __name__ == '__main__':
    main()

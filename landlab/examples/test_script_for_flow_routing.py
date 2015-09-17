#! /usr/env/python
"""
just a little script for testing the d8 flow routing class
and drainage area calculation
"""
from __future__ import print_function

from landlab import RasterModelGrid
from landlab.examples.flowRoutingD8 import RouteFlowD8
from landlab.examples.drainageArea import CalcDrainageArea
from numpy import *


def main():
    nr = 5
    nc = 6
    ncells = nr * nc
    dx = 1
    # instantiate grid
    rg = RasterModelGrid(nr, nc, dx)

    elevations = zeros(ncells)
    # set-up interior elevations with random numbers
    # for i in range(0, ncells):
    #    if rg.is_interior(i):
    #        elevations[i]=random.random_sample()

    # set-up with prescribed elevations to test drainage area calcualtion
    helper = [7, 13, 19, 10, 16, 22]
    for i in range(0, 6):
        # print 'helper[i]', helper[i]
        elevations[helper[i]] = 2
    helper = [8, 14, 20, 21, 15, 9]
    for i in range(0, 6):
        elevations[helper[i]] = 3

    # tried making a pit, drainage area algorithm doesn't crash, but it doesn't
    # work perfectly either.
    # elevations[15]=-50

    # printing elevations for debugging purposes
    print('elevation vector')
    print(elevations)

    # instantiate flow routing variable
    flow = RouteFlowD8(ncells)
    # calculate flow directions
    flow_directions = flow.calc_flowdirs(rg, elevations)
    # printing flow directions for debugging purposes
    print('flow direction vector')
    print(flow_directions)

    # instantiate drainage area variable
    da_calculator = CalcDrainageArea(ncells)
    # calculate drainage area
    drain_area = da_calculator.calc_DA(rg, flow_directions)
    # printing drainage area for debugging purposes
    print('drainage area vector')
    print(drain_area)

    # test of other code bit in raster model grid
    #[ms, ma] = rg.calculate_max_gradient_across_node(elevations,12)
    # print 'max slope', ms
    # print 'max angle', ma


if __name__ == '__main__':
    main()

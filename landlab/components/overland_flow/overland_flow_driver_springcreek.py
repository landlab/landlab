#! /usr/env/python
"""

This tests the overland flow and shear stress generator.

"""

import landlab
from landlab.components.dem_support.dem_boundary_conditions import WatershedBoundaryConditions
from landlab.components.overland_flow.overland_flow_generator import OverlandFlow
from landlab.components.uniform_precip.generate_uniform_precip import PrecipitationDistribution
from landlab.io import read_esri_ascii
#import pylab as pl
import pylab
import numpy as np
import os
import time

def main():
    """
    do some stuff
    """
    
    start_time = time.time()
    
    dem_name = 'ExampleDEM/small_sc.asc'
    
#    outlet_row = 82
#    outlet_column = 38
#    next_to_outlet_row = 81
#    next_to_outlet_column = 38


    # Create and initialize a raster model grid by reading a DEM
    DATA_FILE = os.path.join(os.path.dirname(__file__), dem_name)
    print('Reading data from "'+str(DATA_FILE)+'"')
    (rg, z) = read_esri_ascii(DATA_FILE)
    print('DEM has '+str(rg.nrows)+' rows, '+str(rg.ncols)+ \
          ' columns, and cell size '+str(rg.dx))
    #instantiate grid

#    
#    # Modify the grid DEM to set all nodata nodes to inactive boundaries
#    nodata_val=0
#    rg.deactivate_nodata_nodes(z, nodata_val) # set nodata nodes to inactive bounds
#
#    # Set the open boundary (outlet) cell. We want to remember the ID of the 
#    # outlet node and the ID of the interior node adjacent to it. We'll make
#    # the outlet node an open boundary.
#    outlet_node = rg.grid_coords_to_node_id(outlet_row, outlet_column)
#    node_next_to_outlet = rg.grid_coords_to_node_id(next_to_outlet_row, 
#                                                    next_to_outlet_column)
#    rg.set_fixed_value_boundaries(outlet_node)
#            
#    #elevations  = nodata_val*np.ones( nnodes )    
#    ##set-up interior elevations with random numbers
#    ##for i in range(0, nnodes):
#    ##    if rg.is_interior(i):
#    ##        elevations[i]=random.random_sample()
#    #
#    ##set-up with prescribed elevations to test drainage area calcualtion
#    #helper = [7,8,9,10,13,14,15,16]
#    #elevations[helper]=2
#    #helper = [19,20,21,22]
#    #elevations[helper]=3        
#    #elevations[7]=1    
#    #
#    # Get a 2D array version of the elevations
#    elev_raster = rg.node_vector_to_raster(z,True)
#    
#    # Plot topography
#    #pylab.figure(22)
#    #pylab.subplot(121)
#    #im = pylab.imshow(elev_raster, cmap=pylab.cm.RdBu, extent=[0, rg.ncols*rg.dx, 0, rg.nrows*rg.dx])
#    #cb = pylab.colorbar(im)
#    #cb.set_label('Elevation (m)', fontsize=12)
#    #pylab.title('Topography')
#    #
#    #pylab.show()
#    
#    of=OverlandFlow('input_data.txt',rg,0)
#    #rainfall = PrecipitationDistribution()
#    #rainfall.initialize('input_data.txt')
#    #rainfall.update
#    #
#    ###for now this is in hours, so put into seconds
#    #storm_duration = rainfall.storm_duration*3600
#    ###in mm/hour, so convert to m/second
#    #storm_intensity = rainfall.intensity/1000/3600
#    #print "storm duration, seconds ", storm_duration
#    #print "storm duration, hours ", rainfall.storm_duration
#    #print "storm intensity ", storm_intensity
#    
#    #tau = of.run_one_step(rg,z,outlet_node,node_next_to_outlet,storm_duration,storm_intensity)
#    of.run_one_step(rg,z,outlet_node,node_next_to_outlet,60*40,2.78e-05, 60*15)
#    
#    #imshow(elev_raster)
#    #colorbar()
#    #show()
#    ##
#    ### Create a shaded image
#    ##pylab.close()  # clear any pre-existing plot
#    ##image_extent = [0, 0.001*dx*nc, 0, 0.001*dx*nr] # in km
#    ##im = pylab.imshow(elev_raster, cmap=pylab.cm.RdBu, extent=image_extent)
#    ##pylab.xlabel('Distance (km)', fontsize=12)
#    ##pylab.ylabel('Distance (km)', fontsize=12)
#    
#    print('Total run time = '+str(time.time()-start_time)+' seconds.')
    
if __name__ == "__main__":
    main()
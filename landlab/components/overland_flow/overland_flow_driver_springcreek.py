#! /usr/env/python
"""

This tests the overland flow and shear stress generator.

"""

import landlab
#NG still trying to get the waterhsed boundary conditions working
#from landlab.components.dem_support.dem_boundary_conditions import WatershedBoundaryConditions
from landlab.components.overland_flow.overland_flow_generator import OverlandFlow
from landlab.components.uniform_precip.generate_uniform_precip import PrecipitationDistribution
from landlab.io import read_esri_ascii
import pylab
import numpy as np
import os
import time

def main():
    """
    initialize DEM
    Get storm info
    Generate hydrograph
    """
    
    start_time = time.time()
    
    dem_name = 'ExampleDEM/HalfFork.asc'
    



    # Create and initialize a raster model grid by reading a DEM
    DATA_FILE = os.path.join(os.path.dirname(__file__), dem_name)
    print('Reading data from "'+str(DATA_FILE)+'"')
    (rg, z) = read_esri_ascii(DATA_FILE)
    nodata_val=-9999
    
    print('DEM has ' +
          str(rg.number_of_node_rows) + ' rows, ' +
          str(rg.number_of_node_columns) + ' columns, and cell size ' +
          str(rg.dx))
    
    #Below for finding outlet, but not working. boo.
    #wbc=WatershedBoundaryConditions()
    #wbc.set_bc_find_outlet(rg, z, nodata_val)
      
    #Below points are actually the points of interest for plotting the hydrograph
    #They do not have to be THE outlet of the watershed.
    outlet_row = 238
    outlet_column = 216
    next_to_outlet_row = 239
    next_to_outlet_column = 216
          
    outlet_node = rg.grid_coords_to_node_id(outlet_row, outlet_column)
    node_next_to_outlet = rg.grid_coords_to_node_id(next_to_outlet_row, 
                                                    next_to_outlet_column)
                                                                                                       
    # Get a 2D array version of the elevations for plotting purposes
    elev_raster = rg.node_vector_to_raster(z,True)
    
    # Plot topography
    #pylab.figure(22)
    #pylab.subplot(121)
    #im = pylab.imshow(elev_raster, cmap=pylab.cm.RdBu, extent=[0, rg.number_of_node_columns*rg.dx, 0, rg.number_of_node_rows*rg.dx])
    #cb = pylab.colorbar(im)
    #cb.set_label('Elevation (m)', fontsize=12)
    #pylab.title('Topography')
    #
    #pylab.show()
    
    of=OverlandFlow('input_data.txt',rg,0)
    
    #Below for using rainfall distribution
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
    #of.run_one_step(rg,z,outlet_node,node_next_to_outlet,60*15,1.39e-05, 60*15)
    ##^^^^^ REWRITE WITH SEVERAL NODES TO SPEED UP FOR AGU

#    print('Total run time = '+str(time.time()-start_time)+' seconds.')
    
if __name__ == "__main__":
    main()

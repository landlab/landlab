#! /usr/env/python
"""

This tests the overland flow and shear stress generator.

"""

import landlab
from landlab.components.overland_flow.overland_flow_generator_jma import OverlandFlow
from landlab.components.uniform_precip.generate_uniform_precip import PrecipitationDistribution
from landlab.io import read_esri_ascii
from matplotlib import pyplot as plt
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
    
    dem_name = 'HalfFork.asc'
    input_file = 'input_data.txt'
    IT_FILE = os.path.join(os.path.dirname(__file__), input_file)

    # Create and initialize a raster model grid by reading a DEM
    DATA_FILE = os.path.join(os.path.dirname(__file__), dem_name)
    print('Reading data from "'+str(DATA_FILE)+'"')
    (rg, z) = read_esri_ascii(DATA_FILE)
    nodata_val=-9999
    # Modify the grid DEM to set all nodata nodes to inactive boundaries
    rg.set_nodata_nodes_to_inactive(z, nodata_val) # set nodata nodes to inactive bounds
    
    print('DEM has ' +
          str(rg.number_of_node_rows) + ' rows, ' +
          str(rg.number_of_node_columns) + ' columns, and cell size ' +
          str(rg.dx))
    

    # Select point to sample at.
    
    study_row = 97
    study_column = 164
      
    study_node = rg.grid_coords_to_node_id(study_row, study_column)


    ## Set outlet point to set boundary conditions.
    the_outlet_row = 240
    the_outlet_column = 215
    the_outlet_node = rg.grid_coords_to_node_id(the_outlet_row, the_outlet_column)

    rg.set_fixed_value_boundaries(the_outlet_node)
                                                                                                       
    # Get a 2D array version of the elevations for plotting purposes
    elev_raster = rg.node_vector_to_raster(z,True)
    
    # Everything below plots the topography and sampling points
    
    levels = []
    x_up = 1990
    while x_up !=2200:
        levels.append(x_up)
        x_up+=1
    s = plt.contourf(elev_raster, levels, colors='k')#('r','g','b'))
    plt.set_cmap('bone')
    cb = plt.colorbar()

    plt.plot([164], [152],'rs', label ='NE Tributary')
    plt.plot([210],[111], 'ms', label = 'E Tributary')
    plt.plot([150],[109],'cs', label= 'Upper Main Channel')
    plt.plot([215],[9], 'wo', label= 'Outlet')
    plt.legend(loc=3)


    of=OverlandFlow(IT_FILE,rg,0)

    ## (-,-,-,-,how long we route overland flow, intensity in m/s, duration of storm) ##
    ##    of.run_one_step(rg,z,study_node,60*15,100./1000./3600., 60*15)

    ## Trial 1, 5 year storm ##
    of.run_one_step(rg,z,study_node,7500,(7.88*(10**-6)), 2232)
    ## Trial 1, 10 year storm ##
    #of.run_one_step(rg,z,study_node,9000,(7.167*(10**-6)), 2916)
    ## Trial 1, 50 year storm ##
    #of.run_one_step(rg,z,study_node,17500,(7.88*(10**-6)), 17500)



    plt.show()


if __name__ == "__main__":
    main()

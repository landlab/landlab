#! /usr/env/python
"""

This tests the overland flow and shear stress generator.

"""

import landlab
from landlab.components.overland_flow.generate_overland_flow_DEM import OverlandFlow
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
    
    study_row = 110
    study_column = 150
      
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

    plt.plot([150],[109],'cs', label= 'Study Node')
    plt.plot([215],[9], 'wo', label= 'Outlet')
    plt.legend(loc=3)


    of=OverlandFlow()#IT_FILE,rg,0)
    of.initialize(rg)

    ## (-,-,-,-,how long we route overland flow, intensity in m/s, duration of storm) ##

    ## Trial 1, 10 year storm ##
    of.calculate_flow_at_one_point(rg,z,study_node,90,(7.167*(10**-6)), 2916)


    plt.show()


if __name__ == "__main__":
    main()

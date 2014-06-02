#! /usr/env/python
"""

This tests the overland flow and shear stress generator.

"""

import landlab
from landlab.components.overland_flow.generate_overland_flow_DEM1 import OverlandFlow
#from landlab.components.uniform_precip.generate_uniform_precip import PrecipitationDistribution
from landlab.io import read_esri_ascii
from matplotlib import pyplot as plt
#import numpy as np
import os
import time





def main():
    """
    initialize DEM
    Get storm info
    Generate hydrograph
    """
    
    start_time = time.time()
    
    dem_name = 'chiri.asc'
    #input_file = 'input_data.txt'
    #IT_FILE = os.path.join(os.path.dirname(__file__), input_file)

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
    
    study_row = 122
    study_column = 140
      
    study_node = rg.grid_coords_to_node_id(study_row, study_column)


    ## Set outlet point to set boundary conditions.
    the_outlet_row = 152
    the_outlet_column = 230
    the_outlet_node = rg.grid_coords_to_node_id(the_outlet_row, the_outlet_column)

    rg.set_fixed_value_boundaries(the_outlet_node)
                                                                                                       
    # Get a 2D array version of the elevations for plotting purposes
    elev_raster = rg.node_vector_to_raster(z,True)
    
    # Everything below plots the topography and sampling points
    
#    levels = []
#    x_up = 2475
#    while x_up !=2600:
#        levels.append(x_up)
#        x_up+=1
#    plt.figure('Topography')
#    plt.contourf(elev_raster, levels, colors='k')#('r','g','b'))
#    plt.set_cmap('bone')
#    plt.colorbar()
#
#    plt.plot([219],[85],'cs', label= 'Study Node')
#    plt.plot([224],[75], 'wo', label= 'Outlet')
#    plt.legend(loc=3)
    
    #pd = PrecipitationDistribution()
    #pd.initialize()
    #duration_hrs = pd.storm_duration
    #intensity_mmhr = pd.intensity

    #duration_secs = duration_hrs*60.0*60.0
    #intensity_ms = ((intensity_mmhr/1000.0)/3600.0)
    #total_duration_secs = 1.25 * duration_secs
    
    #print 'total_duration_secs: ', total_duration_secs


    of=OverlandFlow()
    of.initialize(rg)

    #of.flow_at_one_node(rg, z, study_node, 5868, (7.064*(10**-6)), 5868)
    #of.plot_at_one_node()

    of.flow_across_grid(rg, z, 5800, (9.2177*(10**-6)), 5868)
    of.plot_water_depths(rg)
    of.plot_discharge(rg)
    of.plot_shear_stress_grid(rg)
    of.plot_slopes(rg)

    endtime = time.time()
    print endtime - start_time, "seconds"
    plt.show()
    plt.show()


if __name__ == "__main__":
    main()

#! /usr/env/python
""" overland_flow_driver_chiri.py

This is a sample driver which utilizes the
OverlandFlow class from generate_overland_flow_DEM.py
across a subwatershed in Chiricahua Mountains, Arizona.

Written by Jordan Adams, Greg Tucker and Nicole Gasparini.

"""

from landlab.components.overland_flow.generate_overland_flow_DEM import OverlandFlow
from landlab.io import read_esri_ascii
from matplotlib import pyplot as plt
from landlab.plot import imshow_grid
import os
import numpy as np
import time


def main():
    """
    This driver takes in a DEM of a subwatershed from
    the Chiricahua Mountains, Arizona and routes a storm across it using
    the default input file (overland_flow_input.txt).

    This has two ways to look at the data: at one point (generating
    a hydrograph and looking for temporal changes in water depth, etc...)
    and across the grid for spatial patterns. 
    
    """
    
    # This provides us with an initial time. At the end, it gives us total
    # model run time in seconds.
    start_time = time.time()
    
    # This is the DEM of the subwatershed from the Chiricahua Mountains, Arizona
    dem_name = 'chiri_10.asc'
    
    # Now we can create and initialize a raster model grid by reading a DEM
    
    # First, this looks for the DEM in the overland_flow folder in Landlab
    DATA_FILE = os.path.join(os.path.dirname(__file__), dem_name)

    # This print statement verifies that we are opening the data file.
    print('Reading data from "'+str(DATA_FILE)+'"')
    
    # Now the ASCII is read, assuming that it it standard ESRI format.
    (rg, z) = read_esri_ascii(DATA_FILE)
    
    # Whatever the NODATA value is in the DEM is needed here to set bounary condition.
    nodata_val=-9999
    
    # Modify the grid DEM to set all NODATA nodes to inactive boundaries
    rg.set_nodata_nodes_to_inactive(z, nodata_val) 
    
    # This prints standard grid characteristics (rows, columns and cell size)

    print('DEM has ' +
          str(rg.number_of_node_rows) + ' rows, ' +
          str(rg.number_of_node_columns) + ' columns, and cell size ' +
          str(rg.dx))
    


    # Right now, the outlet must be explicitly set for boundary conditions
    # using the row and column from the DEM.
    #CHIRICAHUA
    the_outlet_row = 1
    the_outlet_column = 29
    
    sampling_row = 2
    sampling_col = 29

    
    # This converts the grid row and column into coordinates that the raster grid will recognize
    the_outlet_node = rg.grid_coords_to_node_id(the_outlet_row, the_outlet_column)

    # The outlet node is set as a fixed value boundary
    rg.set_fixed_value_boundaries(the_outlet_node)
                                                                                                       
    # Now we will initialize the Overland Flow component.
    of=OverlandFlow(rg)   

    # First, the flow_at_one_node() method will be covered here.
    
    # When using the flow_at_one_node() method, a study node is needed to sample at.
    # This should not be a boundary node!        

    
    # This takes the study row and study column grid coordinates and converts it 
    # to a node ID that the raster grid will recognize. 
    study_node = rg.grid_coords_to_node_id(sampling_row, sampling_col)

    # Because this function reads in data using the Model Parameter Dictionary
    # from the default input file, the only arguments needed to run the flow_at_one_node()
    # method is the RasterModelGrid instance, the initial elevations and the study node
    # coordinates.
    
    # This is the standard way to call the flow_at_one_node method using the input file.
    
    #of.flow_at_one_node(rg, z)

    # Once run, we can plot the output
    
    #of.plot_at_one_node()
    
    # If you did not want to use the input file and instead define the total model run time of 
    # rainfall intensity and storm duration, the commented function call below demonstrates this with
    # the following parameters:
        # Total model run time: 7500 seconds
        # Storm intensity: (9.2177*(10**-6)) meters per second.
        # Storm duration: 5868 seconds
    #of.flow_at_one_node(rg, z, study_node, rainfall_duration=5868, rainfall_intensity=(9.2177*(10**-6)), model_duration=7500)


    # Now the flow_across_grid() method will be discussed. The commands needed to
    # run this method are triple-commented ('###') for convenience. All flow_at_one_node()
    # methods calls *should* be commented out to run this driver quicker.
    
    ###of.flow_across_grid(rg, z)
    # Again, this function reads in data using the ModelParameterDictionary and the default
    # input file. The only required arguments for this method are the  RasterModelGrid instance
    # and the the initial elevations.

    
    # To forgo the call to the input file and define model run time, rainfall
    # intensity and storm duration in the function call itself, the following command 
    # can be used with the following parameters:
        # Total model run time: 7500 seconds
        # Storm intensity: (9.2177*(10**-6)) meters per second.
        # Storm duration: 5868 seconds
    of.flow_across_grid(rg, z, rainfall_duration=5868, rainfall_intensity=(9.2177*(10**-3)), model_duration=7500)

    endtime = time.time()
    print endtime - start_time, "seconds"
    
    
    ## Plotting topography
    z[np.where(z<=0.)] = 9999            # temporarily change their elevs ...
    zmin = np.amin(z)                    # ... so we can find the minimum ...
    z[np.where(z==9999)] = zmin          # ... and assign them this value.
    
    plt.figure('Topography')
    imshow_grid(rg, z, show=True)
    
    # Plotting water depths
    plt.figure('Water Depths')
    imshow_grid(rg, of.h, show=True)



if __name__ == "__main__":
    main()

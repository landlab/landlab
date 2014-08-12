#! /usr/env/python
""" overland_flow_driver_chiri.py

This is a sample driver which utilizes the
OverlandFlow class from generate_overland_flow_DEM.py
across a subwatershed in Chiricahua Mountains, Arizona.

Written by Jordan Adams, Greg Tucker and Nicole Gasparini.

"""

from landlab.components.overland_flow.generate_overland_flow_DEM_fields import OverlandFlow
from landlab.io import read_esri_ascii
from matplotlib import pyplot as plt
from landlab.plot import imshow_field
import os
import time

def plot_topography(grid, elev):
    ''' 
    This function takes the DEM read in below and plots the
    elevations for visualization purposes.
    '''
    # Get a 2D array version of the elevations for plotting purposes
    elev_raster = grid.node_vector_to_raster(elev,True)
    
    # Everything below plots the topography and sampling points
    levels = []
    # To better create a colorbar...
    x_up = 2475
    while x_up !=2600:
        levels.append(x_up)
        x_up+=1
    plt.figure('Topography')
    plt.contourf(elev_raster, levels, colors='k')
    plt.set_cmap('bone')
    plt.colorbar()
    
    # To plot the study node and outlet node on the DEM...
    plt.plot([122],[140],'cs', label= 'Study Node')
    plt.plot([152],[230], 'wo', label= 'Outlet')
    
    plt.legend(loc=3)

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
    dem_name = 'chiri.asc'
    
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
    the_outlet_row = 152
    the_outlet_column = 230
    
    # This converts the grid row and column into coordinates that the raster grid will recognize
    the_outlet_node = rg.grid_coords_to_node_id(the_outlet_row, the_outlet_column)

    # The outlet node is set as a fixed value boundary
    rg.set_fixed_value_boundaries(the_outlet_node)
                                                                                                       
    # To plot the grid, we can call the plot_topography() function
    plot_topography(rg, z)
    
    # Now we will initialize the Overland Flow component.
    of=OverlandFlow(rg)   

    # First, the flow_at_one_node() method will be covered here.
    
    # When using the flow_at_one_node() method, a study node is needed to sample at.
    # This should not be a boundary node!        

    study_row = 122
    study_column = 140
    
    # This takes the study row and study column grid coordinates and converts it 
    # to a node ID that the raster grid will recognize. 
    study_node = rg.grid_coords_to_node_id(study_row, study_column)

    # Because this function reads in data using the Model Parameter Dictionary
    # from the default input file, the only arguments needed to run the flow_at_one_node()
    # method is the RasterModelGrid instance, the initial elevations and the study node
    # coordinates.
    
    # This is the standard way to call the flow_at_one_node method using the input file.
    of.flow_at_one_node(rg, z, study_node)

    # Now we can update the function, using output from the above function call as the initial condition for our update
    of.update_at_one_point(rg, rainfall_duration=900, model_duration=900,rainfall_intensity=0.00000193333)#
    of.update_at_one_point(rg, rainfall_duration=480, model_duration=480,rainfall_intensity=0.00000183333)# 

    # Once run, we can plot the output
    of.plot_at_one_node()
    of.plot_water_depths(rg)
    

    # If you did not want to use the input file and instead define the total model run time of 
    # rainfall intensity and storm duration, the commented function call below demonstrates this with
    # the following parameters:
        # Total model run time: 6000 seconds
        # Storm intensity: 0.0000137778 meters per second.
        # Storm duration: 900 seconds
    #of.flow_at_one_node(rg, z, study_node, rainfall_duration=900, rainfall_intensity=0.0000137778,model_duration=900)
    #of.update_at_one_point(rg, rainfall_duration=900, model_duration=900,rainfall_intensity=0.00000193333)#
    #of.update_at_one_point(rg, rainfall_duration=480, model_duration=480,rainfall_intensity=0.00000183333)#
    
    
    # Now the flow_across_grid() method will be discussed. The commands needed to
    # run this method are triple-commented ('###') for convenience. All flow_at_one_node()
    # methods calls *should* be commented out to run this driver quicker.
    
    # Again, this function reads in data using the ModelParameterDictionary and the default
    # input file. The only required arguments for this method are the  RasterModelGrid instance
    # and the the initial elevations.
    
    # This is the standard call to the flow_across_grid() method
    # To model this across the grid, uncomment out the lines below...
    ###of.flow_across_grid(rg, z)
    
    # And to update across the grid...
    ###of.update_across_grid(rg, rainfall_duration=900, model_duration=900,rainfall_intensity=0.00000193333)#
    ###of.update_across_grid(rg, rainfall_duration=480, model_duration=480,rainfall_intensity=0.00000183333)#
    
    # And these are the calls to plot the rasters
    of.plot_water_depths(rg)
    of.plot_slopes(rg)

    ###of.plot_shear_stress_grid(rg)
    

    endtime = time.time()
    print endtime - start_time, "seconds"
    
    # And to test our fields capability, we'll print the largest value of each...
    d = rg['node']['water_depth']
    print 'max water depth: ', max(d)
    
    t = rg['node']['shear_stress']
    print 'max shear stress: ', max(t)
    
    q = rg['link']['water_discharge']
    print 'max discharge: ', max(q)
    
    
    #If you call the plot_topography method...
    plt.show()


if __name__ == "__main__":
    main()

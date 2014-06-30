#! /usr/env/python
""" overland_flow_driver_springcreek.py
<<<<<<< HEAD

This is a sample driver which utilizes the
OverlandFlow class from generate_overland_flow_DEM.py
across a subwatershed in Spring Creek, Colorado.
=======

 This component simulates overland flow using
 the 2-D numerical model of shallow-water flow
 over topography using the Bates et al. (2010)
 algorithm for storage-cell inundation modeling
 across a small watershed in Spring Creek, Colorado

Written by Jordan Adams, Nicole Gasparini and Greg Tucker

>>>>>>> origin/master

"""

from landlab.components.overland_flow.generate_overland_flow_DEM import OverlandFlow
from landlab.io import read_esri_ascii
from matplotlib import pyplot as plt
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
    x_up = 1990
    while x_up !=2200:
        levels.append(x_up)
        x_up+=1
    plt.figure('Topography')
    plt.contourf(elev_raster, levels, colors='k')#('r','g','b'))
    plt.set_cmap('bone')
    plt.colorbar()
    
    # To plot the study node and outlet node on the DEM...
    #plt.plot([150],[109],'cs', label= 'Study Node')
    #plt.plot([215],[9], 'wo', label= 'Outlet')
    plt.legend(loc=3)

def main():
    """
    This driver takes in a DEM of a subwatershed from
    Spring Creek, Colorado and routes a storm across it using
    the default input file (overland_flow_input.txt).

    This has two ways to look at the data: at one point (generating
    a hydrograph and looking for temporal changes in water depth, etc...)
    and across the grid for spatial patterns. 
    
    
    """
    # This provides us with an initial time. At the end, it gives us total
    # model run time in seconds.
    start_time = time.time()
    
    # This is the DEM of the subwatershed from Spring Creek, Colorado
    dem_name = 'HalfFork.asc'

    # Now we can create and initialize a raster model grid by reading a DEM
    
    # First, this looks for the DEM in the overland_flow folder in Landlab
    DATA_FILE = os.path.join(os.path.dirname(__file__), dem_name)
    
    # This print statement verifies that we are opening the data file.
    print('Reading data from "'+str(DATA_FILE)+'"')
    
    # Now the ASCII is read, assuming that it is standard ESRI format.
    (rg, z) = read_esri_ascii(DATA_FILE)
    
    # Whatever the NODATA value is in the DEM is needed here to set bounary condition.
    nodata_val=-9999
    
    # Modify the grid DEM to set all NODATA nodes to inactive boundaries
    rg.set_nodata_nodes_to_inactive(z, nodata_val) 
    
    # This gives standard grid characteristics (rows, columns and cell size)
    print('DEM has ' +
          str(rg.number_of_node_rows) + ' rows, ' +
          str(rg.number_of_node_columns) + ' columns, and cell size ' +
          str(rg.dx))


    # Right now, the outlet must be explicitly set for boundary conditions
    # using the row and column from the DEM.
    the_outlet_row = 240
    the_outlet_column = 215
    
    # This converts the grid row and column into coordinates that the raster grid will recognize
    the_outlet_node = rg.grid_coords_to_node_id(the_outlet_row, the_outlet_column)

    # The outlet node is set as a fixed value boundary
    rg.set_fixed_value_boundaries(the_outlet_node)
                                                                                                       
    # To plot the grid, we can call the plot_topography() function
    #plot_topography(rg, z)
    
    # Now we will initialize the Overland Flow component.
    of=OverlandFlow(rg)   
    
    # First, the flow_at_one_node() method will be covered here.
    
    # When using the flow_at_one_node() method, a study node is needed to sample at.
    # This should not be a boundary node!    
    study_row = 110
    study_column = 150
    
    # This takes the study row and study column grid coordinates and converts it 
    # to a node ID that the raster grid will recognize.  
    study_node = rg.grid_coords_to_node_id(study_row, study_column)

    # Because this function reads in data using the Model Parameter Dictionary
    # from the default input file, the only arguments needed to run the flow_at_one_node()
    # method is the RasterModelGrid instance, the initial elevations and the study node
    # coordinates.
    
    # This is the standard way to call the flow_at_one_node method using the input file.
    of.flow_at_one_node(rg, z, study_node, 7500)

    # Once run, we can plot the output
    of.plot_at_one_node()
    
    # If you did not want to use the input file and instead define the total model run time of 
    # rainfall intensity and storm duration, the commented function call below demonstrates this with
    # the following parameters:
        # Total model run time: 6000 seconds
        # Storm intensity: (9.2177*(10**-6)) meters per second.
        # Storm duration: 5868 seconds
    #of.flow_at_one_node(rg, z, study_node, 6000, (9.2177*(10**-6)),5868)

    # Now the flow_across_grid() method will be discussed. The commands needed to
    # run this method are triple-commented ('###') for convenience. All flow_at_one_node()
    # methods calls *should* be commented out to run this driver quicker.
    
    # Again, this function reads in data using the ModelParameterDictionary and the default
    # input file. The only required arguments for this method are the  RasterModelGrid instance
    # and the the initial elevations.
    
    # This is the standard call to the flow_across_grid() method
    ###of.flow_across_grid(rg,z)
    
    # And these are the calls to plot the rasters
    ###of.plot_water_depths(rg)
    ###of.plot_discharge(rg)
    ###of.plot_shear_stress_grid(rg)
    ###of.plot_slopes(rg)
    
    # To forgo the call to the input file and define model run time, rainfall
    # intensity and storm duration in the function call itself, the following command 
    # can be used with the following parameters:
        # Total model run time: 6000 seconds
        # Storm intensity: (9.2177*(10**-6)) meters per second.
        # Storm duration: 5868 seconds
    ###of.flow_across_grid(rg, z, 6000, (9.2177*(10**-6)), 5868)


    endtime = time.time()
    print endtime - start_time, "seconds"
    
    #If you call the plot_topography method...
    #plt.show()


if __name__ == "__main__":
    main()

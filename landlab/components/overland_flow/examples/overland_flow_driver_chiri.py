#! /usr/env/python
""" overland_flow_driver_chiri.py

This is a sample driver which utilizes the
OverlandFlow class from generate_overland_flow_DEM.py
across a subwatershed in Chiricahua Mountains, Arizona.

Written by Jordan Adams, Greg Tucker and Nicole Gasparini.

"""
from __future__ import print_function

from landlab.components.overland_flow.generate_overland_flow_DEM import OverlandFlow
from landlab.io import read_esri_ascii
from matplotlib import pyplot as plt
import os
import time
from landlab.plot import imshow_grid


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
print(('Reading data from "'+str(DATA_FILE)+'"'))

# Now the ASCII is read, assuming that it it standard ESRI format.
(rg, z) = read_esri_ascii(DATA_FILE)

# Whatever the NODATA value is in the DEM is needed here to set bounary condition.
nodata_val=-9999

# Modify the grid DEM to set all NODATA nodes to inactive boundaries
rg.set_nodata_nodes_to_inactive(z, nodata_val) 

# This prints standard grid characteristics (rows, columns and cell size)

print(('DEM has ' +
        str(rg.number_of_node_rows) + ' rows, ' +
        str(rg.number_of_node_columns) + ' columns, and cell size ' +
        str(rg.dx)))

# Right now, the outlet must be explicitly set for boundary conditions
# using the row and column from the DEM.
the_outlet_row = 1
the_outlet_column = 29

# This converts the grid row and column into coordinates that the raster grid will recognize
the_outlet_node = rg.grid_coords_to_node_id(the_outlet_row, the_outlet_column)

# The outlet node is set as a fixed value boundary
rg.set_fixed_value_boundaries(the_outlet_node)
                                                                                                
# Now we will initialize the Overland Flow component.
of=OverlandFlow(rg)   

# First, the flow_at_one_node() method will be covered here.
# When using the flow_at_one_node() method, a study node is needed to sample at.
# This should not be a boundary node!        

# To use the flow_at_one_node() method, uncomment out every line with ('##')

##study_row = 5
##study_column = 26

# This takes the study row and study column grid coordinates and converts it 
# to a node ID that the raster grid will recognize. 

##study_node = rg.grid_coords_to_node_id(study_row, study_column)

# Because this function reads in data using the Model Parameter Dictionary
# from the default input file, the only arguments needed to run the flow_at_one_node()
# method is the RasterModelGrid instance, the initial elevations and the study node
# coordinates.

# This is the standard way to call the flow_at_one_node method using the input file.

##of.flow_at_one_node(rg, z, study_node)
##of.update_at_one_point(rg, rainfall_duration=900, model_duration=900,rainfall_intensity=0.0000193333)
##of.update_at_one_point(rg, rainfall_duration=480, model_duration=480,rainfall_intensity=0.0000183333)

# Once run, we can plot the output

##of.plot_at_one_node()

# Now the flow_across_grid() method will be discussed. The commands needed to
# run this method are triple-commented ('###') for convenience. All flow_at_one_node()
# methods calls *should* be commented out to run this driver quicker.

# Again, this function reads in data using the ModelParameterDictionary and the default
# input file. The only required arguments for this method are the  RasterModelGrid instance
# and the the initial elevations.

# This is the standard call to the flow_across_grid() method
of.flow_across_grid(rg,z)
of.update_across_grid(rg, rainfall_duration=900, model_duration=900,rainfall_intensity=0.0000193333)
of.update_across_grid(rg, rainfall_duration=480, model_duration=900,rainfall_intensity=0.0000183333)

plt.figure('Total Erosion, m')
imshow_grid(rg, of.total_dzdt, show=True)



endtime = time.time()
print(endtime - start_time, "seconds")


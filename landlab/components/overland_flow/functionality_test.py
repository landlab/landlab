## Test script to demonstrate resultant errors in raster_funcs.py
##
##
## May 21, 2014
## Jordan Adams

import landlab
from landlab.io import read_esri_ascii
import os
    
    
dem_name = 'HalfFork.asc'
DATA_FILE = os.path.join(os.path.dirname(__file__), dem_name)


print('Reading data from "'+str(DATA_FILE)+'"')
(rg, z) = read_esri_ascii(DATA_FILE)
nodata_val=-9999

rg.set_nodata_nodes_to_inactive(z, nodata_val) # set nodata nodes to inactive bounds

print('DEM has ' +
    str(rg.number_of_node_rows) + ' rows, ' +
    str(rg.number_of_node_columns) + ' columns, and cell size ' +
    str(rg.dx))
    
the_outlet_row = 240
the_outlet_column = 215
the_outlet_node = rg.grid_coords_to_node_id(the_outlet_row, the_outlet_column)

rg.set_fixed_value_boundaries(the_outlet_node)


interior_nodes = rg.get_active_cell_node_ids()
print interior_nodes
elev_at_interior_nodes = []
for i in range(len(interior_nodes)):
    elev_at_interior_nodes.append(z[interior_nodes[i]])
    
## MULTIPLE ERRORS ##

##################################################################################
## ERROR 1: TYPE ERROR: unhashable type: 'numpy.ndarray'.
## Dan suggests this is because I do not need to pass the Raster grid argument ("rg").

#rg.calculate_steepest_descent_across_cell_faces(rg, z)
#rg.calculate_steepest_descent_across_adjacent_cells(rg,z)

##################################################################################
## ERROR 2: INDEX ERROR: index 9223372036854775807 is out of bounds
## I'm not exactly sure WHAT is happening here. My array of elevations ("z") is of length 63488
## This is what happens when I don't pass the grid as an argument.

#rg.calculate_steepest_descent_across_cell_faces(z)
#rg.calculate_steepest_descent_across_adjacent_cells(z)


##################################################################################
## ERROR 3: ASSERTION ERROR: len(args < 2)
## This is the oldest error I was getting. when I pass the active cell node ids array to the function,
## the assertion error is that the number of arguments in this OPTIONAL argument is too small (less than 2). 
## NOTE I ONLY GET THIS ERROR WHEN THE GRID IS ALSO PASSED.

#rg.calculate_steepest_descent_across_cell_faces(rg, z, interior_nodes)
#rg.calculate_steepest_descent_across_adjacent_cells(rg,z, interior_nodes)


## WHEN I DON'T PASS THE GRID I GET ERROR #2, THE INDEX ERROR.

#rg.calculate_steepest_descent_across_cell_faces(z, interior_nodes)
#rg.calculate_steepest_descent_across_adjacent_cells(z, interior_nodes)




#! /usr/env/python

"""
test_script_for_route_flow_dn.py:

Tests and illustrates use of route_flow_dn component.
"""

from landlab.components.flow_routing import FlowRouter
from landlab.io import read_esri_ascii
from landlab.plot.imshow import imshow_node_grid
import os
import pylab

dem_name = '../../../io/tests/data/west_bijou_gully.asc'
outlet_row = 82
outlet_column = 38

# Read in a DEM and set its boundaries
DATA_FILE = os.path.join(os.path.dirname(__file__), dem_name)
(grid, z) = read_esri_ascii(DATA_FILE)
grid.set_nodata_nodes_to_inactive(z, 0) # set nodata nodes to inactive bounds
outlet_node = grid.grid_coords_to_node_id(outlet_row, outlet_column)

# Route flow
flow_router = FlowRouter(grid)
flow_router.route_flow()

# Get a 2D array version of the elevations
ar = grid.node_vector_to_raster(grid['node']['drainage_area'])

numcols = grid.number_of_node_columns
numrows = grid.number_of_node_rows
dx = grid.dx

# Create a shaded image
pylab.close()  # clear any pre-existing plot
im = pylab.imshow(ar, cmap=pylab.cm.RdBu, extent=[0,numcols*dx,0,numrows*dx])
# add contour lines with labels
cset = pylab.contour(ar, extent=[0,numcols*dx,numrows*dx,0])
pylab.clabel(cset, inline=True, fmt='%1.1f', fontsize=10)

# add a color bar on the side
cb = pylab.colorbar(im)
cb.set_label('Drainage area, sq meters')

# add a title and axis labels
pylab.title('DEM')
pylab.xlabel('Distance (m)')
pylab.ylabel('Distance (m)')

# Display the plot
pylab.show()



from landlab import RasterModelGrid
from landlab import ModelParameterDictionary
from landlab.components.sed_trp_shallow_flow.transport_sed_in_shallow_flow import SurfaceFlowTransport

import time
import pylab
import numpy as np

#get the needed properties to build the grid:
inputs = ModelParameterDictionary('./stof_params.txt')
nrows = inputs.read_int('nrows')
ncols = inputs.read_int('ncols')
dx = inputs.read_float('dx')
h_init = inputs.read_float('h_init')
h_boundary = inputs.read_float('h_boundary')
drop_ht = inputs.read_float('drop_ht')
initial_slope = inputs.read_float('initial_slope')
time_to_run = inputs.read_int('run_time')

left_middle_node = ncols*(nrows/2)
z0 = drop_ht+dx*(ncols-1)*initial_slope

mg = RasterModelGrid(nrows, ncols, dx)
mg.set_inactive_boundaries(True, False, True, False)

#create the fields in the grid
mg.create_node_array_zeros('planet_surface__elevation')
mg.create_node_array_zeros('planet_surface__water_depth')

#set the initial water depths
h = mg.create_node_array_zeros() + h_init
h[left_middle_node] = h_boundary
h[left_middle_node-ncols] = h_boundary
h[left_middle_node+ncols] = h_boundary
mg['node'][ 'planet_surface__water_depth'] = h

# Set initial topography
x = mg.get_node_x_coords()
y = mg.get_node_y_coords()
zinit = mg.create_node_array_zeros()
zinit = z0-initial_slope*x
rightside = mg.right_edge_node_ids()
zinit[rightside] = 0.
mg['node']['planet_surface__elevation'] = zinit

# Display a message
print( 'Running ...' )
start_time = time.time()

#instantiate the component:
transport_component = SurfaceFlowTransport(mg, './stof_params.txt')

#perform the loop:
elapsed_time = 0. #total time in simulation
while elapsed_time < time_to_run:
    timestep = transport_component.set_and_return_dynamic_timestep()
    mg = transport_component.transport_sed(elapsed_time)
    elapsed_time += timestep

#Finalize and plot
zm = mg.at_node['planet_surface__elevation']
h = mg.at_node['planet_surface__water_depth']
ddz=zm-zinit
print ddz[np.where(ddz!=0.)]
print np.amax(ddz)
    
# Get a 2D array version of the water depths and elevations
hr = mg.node_vector_to_raster(h)
zr = mg.node_vector_to_raster(zm)
dzr = mg.node_vector_to_raster(ddz)

# Clear previous plots
pylab.figure(1)
pylab.close()
pylab.figure(2)
pylab.close()

# Plot topography
pylab.figure(1)
pylab.subplot(131)
im = pylab.imshow(zr, cmap=pylab.cm.RdBu)  # display a colored image
pylab.colorbar(im)
pylab.title('Topography')
    
# Plot change in topo
pylab.figure(1)
pylab.subplot(132)
im = pylab.imshow(dzr, cmap=pylab.cm.RdBu)  # display a colored image
pylab.colorbar(im)
pylab.title('Topo change')
    
# Plot water depth
pylab.subplot(133)
im2 = pylab.imshow(hr, cmap=pylab.cm.RdBu)  # display a colored image
#pylab.clim(0, 0.25)
pylab.colorbar(im2)
pylab.title('Water depth')
    
# Display the plots
pylab.show()
print('Done.')
print('Total run time = '+str(time.time()-start_time)+' seconds.')

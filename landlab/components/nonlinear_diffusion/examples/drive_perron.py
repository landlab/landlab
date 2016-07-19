from __future__ import print_function

import numpy
from landlab import RasterModelGrid, CLOSED_BOUNDARY
from landlab import ModelParameterDictionary
from landlab.components.nonlinear_diffusion.Perron_nl_diffuse import PerronNLDiffuse
import pylab
import time


inputs = ModelParameterDictionary('./drive_perron_params.txt')
nrows = inputs.read_int('nrows')
ncols = inputs.read_int('ncols')
dx = inputs.read_float('dx')
dt = inputs.read_float('dt')
time_to_run = inputs.read_float('run_time')
# nt needs defining
uplift = inputs.read_float('uplift_rate')
init_elev = inputs.read_float('init_elev')

mg = RasterModelGrid(nrows, ncols, dx)
#mg.set_inactive_boundaries(False, False, False, False)
# mg.set_inactive_boundaries(True,True,True,True)
mg.set_looped_boundaries(True, True)

#create the fields in the grid
mg.add_zeros('topographic__elevation', at='node')
z = mg.zeros(at='node') + init_elev
mg['node'][ 'topographic__elevation'] = z + numpy.random.rand(len(z))/1000.

# Now add a step to diffuse out:
# mg.at_node['topographic__elevation'][mg.active_nodes[:(mg.active_nodes.shape[0]//2.)]]
# += 0.05 #half block uplift

# pylab.figure(1)
# pylab.close()
#elev = mg['node']['topographic__elevation']
#elev_r = mg.node_vector_to_raster(elev)
# pylab.figure(1)
#im = pylab.imshow(elev_r, cmap=pylab.cm.RdBu)
# pylab.show()

# Display a message
print('Running ...')
start_time = time.time()

# instantiate the component:
diffusion_component = PerronNLDiffuse(mg, './drive_perron_params.txt')

# perform the loop:
elapsed_time = 0.  # total time in simulation
while elapsed_time < time_to_run:
    print(elapsed_time)
    if elapsed_time + dt < time_to_run:
        diffusion_component.input_timestep(dt)
    mg.at_node['topographic__elevation'][mg.core_nodes] += uplift * dt
    # mg.at_node['topographic__elevation'][mg.active_nodes[:(mg.active_nodes.shape[0]//2.)]] += uplift*dt #half block uplift
    # mg.at_node['topographic__elevation'][mg.active_nodes] += (numpy.arange(len(mg.active_nodes))) #nodes are tagged with their ID
    # pylab.figure(1)
    # pylab.close()
    #elev = mg['node']['topographic__elevation']
    #elev_r = mg.node_vector_to_raster(elev)
    # pylab.figure(1)
    #im = pylab.imshow(elev_r, cmap=pylab.cm.RdBu)
    # pylab.show()

    mg = diffusion_component.diffuse(mg, elapsed_time)
    elapsed_time += dt

#Finalize and plot
elev = mg['node']['topographic__elevation']
elev_r = mg.node_vector_to_raster(elev)

# Clear previous plots
pylab.figure(1)
pylab.close()

# Plot topography
pylab.figure(1)
im = pylab.imshow(elev_r, cmap=pylab.cm.RdBu)  # display a colored image
print(elev_r)
pylab.colorbar(im)
pylab.title('Topography')

pylab.figure(2)
# display a colored image
im = pylab.plot(dx * numpy.arange(nrows), elev_r[:, int(ncols // 2)])
pylab.title('Vertical cross section')

pylab.show()

print('Done.')
print(('Total run time = ' + str(time.time() - start_time) + ' seconds.'))

from __future__ import print_function

import numpy
from landlab import RasterModelGrid, CLOSED_BOUNDARY
from landlab.components.nonlinear_diffusion import NonlinearDiffuser
import pylab
import time

mg = RasterModelGrid((100,100), 10.)
mg.set_inactive_boundaries(True, False, True, False)


#create the fields in the grid
mg.add_zeros('topographic__elevation', at='node')
z = mg.zeros(at='node') + 0.5
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
nld = NonlinearDiffuser(mg, nonlinear_diffusivity=0.1, critical_slope=0.3)

# perform the loop:
elapsed_time = 0.  # total time in simulation
time_to_run = 1000000.
uplift = 0.0001
dt = 100
while elapsed_time < time_to_run:
    print(elapsed_time)
    # if elapsed_time + dt < time_to_run:
    #     diffusion_component.input_timestep(dt)
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

    nld.run_one_step(dt)
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
nrows = mg.number_of_node_rows
ncols = mg.number_of_node_columns
im = pylab.plot(10. * numpy.arange(nrows), elev_r[:, int(ncols // 2)])
pylab.title('Vertical cross section')

pylab.show()

print('Done.')
print(('Total run time = ' + str(time.time() - start_time) + ' seconds.'))

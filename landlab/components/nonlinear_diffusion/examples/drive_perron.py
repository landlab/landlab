import numpy
from landlab import RasterModelGrid
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
#nt needs defining
uplift = inputs.read_float('uplift')
init_elev = inputs.read_float('init_elev')

mg = RasterModelGrid(nrows, ncols, dx)
mg.set_inactive_boundaries(False, False, False, False)

#create the fields in the grid
mg.create_node_array_zeros('planet_surface__elevation')
z = mg.create_node_array_zeros() + init_elev
mg['node'][ 'planet_surface__elevation'] = z + numpy.random.rand(len(z))/1000.

#pylab.figure(1)
#pylab.close()
#elev = mg['node']['planet_surface__elevation']
#elev_r = mg.node_vector_to_raster(elev)
#pylab.figure(1)
#im = pylab.imshow(elev_r, cmap=pylab.cm.RdBu)
#pylab.show()

# Display a message
print( 'Running ...' )
start_time = time.time()

#Need to test out if the routine applies its own uplift. We'll see... So don't uplift for now.

#instantiate the component:
diffusion_component = PerronNLDiffuse(mg, './drive_perron_params.txt')

#perform the loop:
elapsed_time = 0. #total time in simulation
while elapsed_time < time_to_run:
    print elapsed_time
    if elapsed_time+dt<time_to_run:
        diffusion_component.gear_timestep(dt)
    else:
        diffusion_component.gear_timestep(time_to_run-elapsed_time)
    mg.at_node['planet_surface__elevation'][mg.active_nodes[:(mg.active_nodes.shape[0]//2.)]] += uplift*dt
    
    #pylab.figure(1)
    #pylab.close()
    #elev = mg['node']['planet_surface__elevation']
    #elev_r = mg.node_vector_to_raster(elev)
    #pylab.figure(1)
    #im = pylab.imshow(elev_r, cmap=pylab.cm.RdBu)
    #pylab.show()
    
    mg = diffusion_component.diffuse(elapsed_time)
    elapsed_time += dt

#Finalize and plot
elev = mg['node']['planet_surface__elevation']
elev_r = mg.node_vector_to_raster(elev)

# Clear previous plots
pylab.figure(1)
pylab.close()

# Plot topography
pylab.figure(1)
im = pylab.imshow(elev_r, cmap=pylab.cm.RdBu)  # display a colored image
print elev_r
pylab.colorbar(im)
pylab.title('Topography')

pylab.show()

print('Done.')
print('Total run time = '+str(time.time()-start_time)+' seconds.')

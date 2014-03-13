from landlab.components.flow_routing.route_flow_dn import FlowRouter
from landlab.components.stream_power.fastscape_stream_power import SPEroder
from landlab.components.nonlinear_diffusion.Perron_nl_diffuse import PerronNLDiffuse
from landlab import ModelParameterDictionary

from landlab import RasterModelGrid
import numpy as np
import pylab

#get the needed properties to build the grid:
input_file = './germany_test_params.txt'
inputs = ModelParameterDictionary(input_file)
nrows = inputs.read_int('nrows')
ncols = inputs.read_int('ncols')
dx = inputs.read_float('dx')
leftmost_elev = inputs.read_float('leftmost_elevation')
initial_slope = inputs.read_float('initial_slope')
uplift_rate = inputs.read_float('uplift_rate')
dt = inputs.read_float('dt')
runtime = inputs.read_float('total_time')
nt = int(runtime//dt)
uplift_per_step = uplift_rate * dt

mg = RasterModelGrid(nrows, ncols, dx)
mg.set_inactive_boundaries(False, False, False, False)

#create the elevation field in the grid
mg.create_node_array_zeros('planet_surface__elevation')
z = mg.create_node_array_zeros() + leftmost_elev
z += initial_slope*np.amax(mg.node_y) - initial_slope*mg.node_y
mg['node'][ 'planet_surface__elevation'] = z + np.random.rand(len(z))/100000.

# Display a message
print 'Running ...' 

#instantiate the component:
fr = FlowRouter(mg)
sp = SPEroder(mg, input_file)
diffuse = PerronNLDiffuse(mg, input_file)

#perform the loops:
for i in xrange(nt):
        mg['node']['planet_surface__elevation'] += uplift_per_step
        mg = fr.route_flow(grid=mg)
        mg = sp.erode(mg)
        mg = diffuse.diffuse(mg, i*dt)
        print 'Completed loop ', i

#Finalize and plot
#elev = mg['node']['planet_surface__elevation']
elev = fr.node_water_discharge
elev_r = mg.node_vector_to_raster(elev)
# Clear previous plots
pylab.figure(1)
pylab.close()
pylab.figure(1)
im = pylab.imshow(elev_r, cmap=pylab.cm.RdBu)  # display a colored image
pylab.colorbar(im)
pylab.title('Topography')

elev = mg['node']['planet_surface__elevation']
#elev = fr.node_water_discharge
elev_r = mg.node_vector_to_raster(elev)
# Clear previous plots
pylab.figure(2)
im = pylab.imshow(elev_r, cmap=pylab.cm.RdBu)  # display a colored image
pylab.colorbar(im)
pylab.title('Topography')

print('Done.')

pylab.show()
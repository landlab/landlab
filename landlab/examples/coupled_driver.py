from __future__ import print_function

from landlab.components.flow_routing.route_flow_dn import FlowRouter
from landlab.components.stream_power.fastscape_stream_power import FastscapeEroder
from landlab.components.nonlinear_diffusion.Perron_nl_diffuse import PerronNLDiffuse
from landlab.components.diffusion.diffusion import LinearDiffuser
from landlab import ModelParameterDictionary
from landlab.plot import channel_profile as prf
from landlab.plot.imshow import imshow_node_grid

from landlab import RasterModelGrid
import numpy as np
import pylab

# get the needed properties to build the grid:
input_file = './coupled_params.txt'
inputs = ModelParameterDictionary(input_file)
nrows = inputs.read_int('nrows')
ncols = inputs.read_int('ncols')
dx = inputs.read_float('dx')
leftmost_elev = inputs.read_float('leftmost_elevation')
initial_slope = inputs.read_float('initial_slope')
uplift_rate = inputs.read_float('uplift_rate')
runtime = inputs.read_float('total_time')
dt = inputs.read_float('dt')
nt = int(runtime // dt)
uplift_per_step = uplift_rate * dt

# instantiate the grid object
mg = RasterModelGrid(nrows, ncols, dx)

# create the elevation field in the grid:
# create the field
mg.create_node_array_zeros('topographic__elevation')
z = mg.create_node_array_zeros() + leftmost_elev
z += initial_slope * np.amax(mg.node_y) - initial_slope * mg.node_y
# put these values plus roughness into that field
mg.at_node['topographic__elevation'] = z + np.random.rand(len(z)) / 100000.

# set up grid's boundary conditions (bottom, right, top, left is inactive)
mg.set_closed_boundaries_at_grid_edges(False, True, False, True)

# Display a message
print('Running ...')

# instantiate the components:
fr = FlowRouter(mg)
sp = FastscapeEroder(mg, input_file)
diffuse = PerronNLDiffuse(mg, input_file)
lin_diffuse = LinearDiffuser(grid=mg, input_stream=input_file)

# perform the loops:
for i in xrange(nt):
    # note the input arguments here are not totally standardized between modules
    #mg = diffuse.diffuse(mg, i*dt)
    mg = lin_diffuse.diffuse(dt)
    mg = fr.route_flow()
    mg = sp.erode(mg)
    mg.at_node['topographic__elevation'][mg.core_nodes] += uplift_per_step

    # plot long profiles along channels
    pylab.figure(6)
    profile_IDs = prf.channel_nodes(mg, mg.at_node['topographic__steepest_slope'],
                                    mg.at_node['drainage_area'], mg.at_node['flow_receiver'])
    dists_upstr = prf.get_distances_upstream(mg, len(mg.at_node['topographic__steepest_slope']),
                                             profile_IDs, mg.at_node['links_to_flow_receiver'])
    prf.plot_profiles(dists_upstr, profile_IDs, mg.at_node[
                      'topographic__elevation'])

    print('Completed loop ', i)

print('Completed the simulation. Plotting...')


#Finalize and plot
# Clear previous plots
pylab.figure(1)
pylab.close()
pylab.figure(1)

# display a colored image
imshow_node_grid(mg, 'water__volume_flux', cmap='PuBu')

pylab.figure(2)
imshow_node_grid(mg, 'topographic__elevation')  # display a colored image

elev = mg['node']['topographic__elevation']
elev_r = mg.node_vector_to_raster(elev)
pylab.figure(3)
pylab.plot(mg.dx * np.arange(nrows), elev_r[:, int(ncols // 2)])
pylab.title('N-S cross_section')

pylab.figure(4)
pylab.plot(mg.dx * np.arange(ncols), elev_r[int(nrows // 4), :])
pylab.title('E-W cross_section')

drainage_areas = mg['node']['drainage_area'][mg.core_nodes]
steepest_slopes = mg['node']['topographic__steepest_slope'][mg.core_nodes]
pylab.figure(5)
pylab.loglog(drainage_areas, steepest_slopes, 'x')
pylab.xlabel('Upstream drainage area, m^2')
pylab.ylabel('Maximum slope')

print('Done.')

pylab.show()

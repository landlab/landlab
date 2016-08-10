from __future__ import print_function

from six.moves import range

from landlab.components.flow_routing import FlowRouter
from landlab.components.stream_power import FastscapeEroder
from landlab import ModelParameterDictionary
from landlab.plot import channel_profile as prf

from landlab import RasterModelGrid, FIXED_VALUE_BOUNDARY
import numpy as np
import pylab

from time import time

import random


#get the needed properties to build the grid:
input_file = './drive_sp_params.txt'
inputs = ModelParameterDictionary(input_file)
nrows = inputs.read_int('nrows')
ncols = inputs.read_int('ncols')
dx = inputs.read_float('dx')
uplift_rate = inputs.read_float('uplift_rate')
runtime = inputs.read_float('run_time')
dt = inputs.read_float('dt')
nt = int(runtime//dt)
uplift_per_step = uplift_rate * dt

#instantiate the grid object
mg = RasterModelGrid(nrows, ncols, dx)

boundary_node_list = mg.boundary_nodes

#set up its boundary conditions (bottom, right, top, left is inactive)
mg.set_inactive_boundaries(True, True, True, True)
mg.set_closed_boundaries_at_grid_edges(False, False, False, False)
mg.status_at_node[20] = FIXED_VALUE_BOUNDARY

##create the elevation field in the grid:
#create the field
mg.add_zeros('topographic__elevation', at='node')
z = mg.zeros(at='node') + init_elev
#put these values plus roughness into that field
mg['node'][ 'topographic__elevation'] = z + np.random.rand(len(z))/100000.

# Display a message
print('Running ...')

# MN: Loop over several changes in the outlet position
for t in range(5):

    # Set a new random outlet position
    mg.set_inactive_boundaries(True, True, True, True)
    #mg.set_closed_boundaries_at_grid_edges(True,True,True,True)
    random_boundary_node = random.choice(boundary_node_list)
    mg.status_at_node[random_boundary_node] = FIXED_VALUE_BOUNDARY

    # MN: Set the elevation of that random outlet boundary node to zero
    #mg['node'][ 'topographic__elevation'][random_boundary_node] = 0

    print('Random boundary node',  random_boundary_node)


    #instantiate the components:
    fr = FlowRouter(mg)
    sp = FastscapeEroder(mg, input_file)

    time_on = time()

    #perform the inner time loops:
    for i in range(nt):
        mg['node']['topographic__elevation'][mg.core_nodes] += uplift_per_step
        mg = fr.route_flow()
        mg = sp.erode(mg)

        #plot long profiles along channels
        pylab.figure(6)
        profile_IDs = prf.channel_nodes(mg, mg.at_node['topographic__steepest_slope'],
                mg.at_node['drainage_area'], mg.at_node['flow__upstream_node_order'],
                mg.at_node['flow__receiver_node'])
        dists_upstr = prf.get_distances_upstream(mg, len(mg.at_node['topographic__steepest_slope']),
                profile_IDs, mg.at_node['flow__link_to_receiver_node'])
        prf.plot_profiles(dists_upstr, profile_IDs, mg.at_node['topographic__elevation'])
        # print 'Completed loop ', i



    print('Completed the simulation. Plotting...')

    #Finalize and plot
    elev = fr.node_water_discharge
    elev_r = mg.node_vector_to_raster(elev)
    # Clear previous plots
    #pylab.figure(1)
    #pylab.close()
    #pylab.figure(1)
    #im = pylab.imshow(elev_r, cmap=pylab.cm.RdBu)  # display a colored image
    #pylab.colorbar(im)
    #pylab.title('Water discharge')

    elev = mg['node']['topographic__elevation']
    elev_r = mg.node_vector_to_raster(elev)
    pylab.figure(t)
    im = pylab.imshow(elev_r, cmap=pylab.cm.RdBu)  # display a colored image
    pylab.colorbar(im)
    pylab.title('Topography')

    #pylab.figure(3)
    #im = pylab.plot(mg.dx*np.arange(nrows), elev_r[:,int(ncols//2)])
    #pylab.title('N-S cross_section')
    #
    #pylab.figure(4)
    #im = pylab.plot(mg.dx*np.arange(ncols), elev_r[int(nrows//4),:])
    #pylab.title('E-W cross_section')
    #
    #drainage_areas = mg['node']['drainage_area'][mg.core_nodes]
    #steepest_slopes = mg['node']['steepest_slope'][mg.core_nodes]
    #pylab.figure(5)
    #pylab.loglog(drainage_areas, steepest_slopes, 'x')
    #pylab.xlabel('Upstream drainage area, m^2')
    #pylab.ylabel('Maximum slope')

pylab.show()
time_off = time()

print('Done.')
print('Time: ', time_off-time_on)


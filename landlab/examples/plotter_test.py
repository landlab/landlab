from __future__ import print_function

from landlab.components.flow_routing.route_flow_dn import FlowRouter
from landlab.components.stream_power.fastscape_stream_power import FastscapeEroder
from landlab import ModelParameterDictionary
from landlab.plot import imshow
from landlab.plot.video_out import VideoPlotter

from landlab import RasterModelGrid
import numpy as np
import pylab

from time import time


# get the needed properties to build the grid:
input_file = './stream_power_params.txt'
inputs = ModelParameterDictionary(input_file)
nrows = inputs.read_int('nrows')
ncols = inputs.read_int('ncols')
dx = inputs.read_float('dx')
leftmost_elev = inputs.read_float('leftmost_elevation')
initial_slope = inputs.read_float('initial_slope')
uplift_rate = inputs.read_float('uplift_rate')

runtime = 40.
dt = 1.

nt = int(runtime // dt)
uplift_per_step = uplift_rate * dt

# instantiate the grid object
mg = RasterModelGrid(nrows, ncols, dx)
# set up its boundary conditions (bottom, right, top, left is inactive)
mg.set_inactive_boundaries(False, True, False, True)

# create the elevation field in the grid:
# create the field
mg.create_node_array_zeros('topographic__elevation')
z = mg.create_node_array_zeros() + leftmost_elev
z += initial_slope * np.amax(mg.node_y) - initial_slope * mg.node_y
# put these values plus roughness into that field
mg['node']['topographic__elevation'] = z + np.random.rand(len(z)) / 100000.

# Display a message
print('Running ...')

# instantiate the components:
fr = FlowRouter(mg)
sp = FastscapeEroder(mg, input_file)
vid = VideoPlotter(mg, data_centering='node')

time_on = time()
# perform the loops:
for i in xrange(nt):
    print('loop ', i)
    mg['node']['topographic__elevation'][mg.core_nodes] += uplift_per_step
    mg = fr.route_flow()
    mg = sp.erode(mg)
    #vid.add_frame(mg, 'topographic__elevation')
    vid.add_frame(mg, mg.hillshade(alt=15.), cmap='gray')


print('Completed the simulation. Plotting...')

time_off = time()

#Finalize and plot

elev = mg['node']['topographic__elevation']
#imshow.imshow_node_grid(mg, elev)

print('Done.')
print('Time: ', time_off - time_on)

# pylab.show()

vid.produce_video(override_min_max=(0, 1))
# vid.produce_video()

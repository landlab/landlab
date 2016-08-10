'''
simple_sp_driver.py

A simple driver implementing Braun-Willett flow routing and then a
(non-fastscape) stream power component.
DEJH, 09/15/14
'''
from __future__ import print_function

from landlab.components.flow_routing import (FlowRouter,
                                             DepressionFinderAndRouter)
from landlab.components.stream_power import StreamPowerEroder, FastscapeEroder

import numpy
from landlab import RasterModelGrid
from landlab import ModelParameterDictionary
from landlab.plot.imshow import imshow_node_grid
import pylab
import time

inputs = ModelParameterDictionary('./drive_sp_params.txt')
nrows = inputs.read_int('nrows')
ncols = inputs.read_int('ncols')
dx = inputs.read_float('dx')
dt = inputs.read_float('dt')
time_to_run = inputs.read_float('run_time')
#nt needs defining
uplift = 10.*inputs.read_float('uplift_rate')
init_elev = inputs.read_float('init_elev')

mg = RasterModelGrid(nrows, ncols, dx)

#create the fields in the grid
mg.add_zeros('topographic__elevation', at='node')
z = mg.zeros(at='node') + init_elev
#z += mg.node_x*0.001
mg['node'][ 'topographic__elevation'] = z + numpy.random.rand(len(z))/1000.

#make some K values in a field to test
mg.at_node['K_values'] = 0.1+numpy.random.rand(nrows*ncols)/10.

print( 'Running ...' )
time_on = time.time()

#instantiate the components:
fr = FlowRouter(mg)
sp = StreamPowerEroder(mg, './drive_sp_params.txt')
lf = DepressionFinderAndRouter(mg)
#load the Fastscape module too, to allow direct comparison
fsp = FastscapeEroder(mg, './drive_sp_params.txt')

#perform the loop:
elapsed_time = 0. #total time in simulation
while elapsed_time < time_to_run:
#for i in range(10):
    print(elapsed_time)
    if elapsed_time+dt>time_to_run:
        print("Short step!")
        dt = time_to_run - elapsed_time
    mg = fr.route_flow(method='D8')
    lf.map_depressions()
    #print 'Area: ', numpy.max(mg.at_node['drainage_area'])
    #mg = fsp.erode(mg)
    mg,_,_ = sp.erode(mg, dt, node_drainage_areas='drainage_area', slopes_at_nodes='topographic__steepest_slope', K_if_used='K_values')
    #add uplift
    mg.at_node['topographic__elevation'][mg.core_nodes] += uplift*dt
    elapsed_time += dt

time_off = time.time()
print('Elapsed time: ', time_off-time_on)

#Finalize and plot
elev = mg['node']['topographic__elevation']
elev_r = mg.node_vector_to_raster(elev)

# Clear previous plots
pylab.figure(1)
pylab.close()

# Plot topography
pylab.figure(1)
im = imshow_node_grid(mg, 'topographic__elevation')  # display a colored image
print(elev_r)

pylab.figure(2)
im = pylab.plot(dx*numpy.arange(nrows), elev_r[:,int(ncols//2)])  # display a colored image
pylab.title('Vertical cross section')

pylab.show()

print('Done.')

'''
simple_sp_driver.py

A simple driver implementing Braun-Willett flow routing and then a 
(non-fastscape) stream power component.
DEJH, 09/15/14
'''

from landlab.components.flow_routing.route_flow_dn import FlowRouter
from landlab.components.stream_power.stream_power import StreamPowerEroder

import numpy
from landlab import RasterModelGrid
from landlab import ModelParameterDictionary
import pylab
import time

inputs = ModelParameterDictionary('./drive_sp_params.txt')
nrows = inputs.read_int('nrows')
ncols = inputs.read_int('ncols')
dx = inputs.read_float('dx')
dt = inputs.read_float('dt')
time_to_run = inputs.read_float('run_time')
#nt needs defining
uplift = inputs.read_float('uplift_rate')
init_elev = inputs.read_float('init_elev')

mg = RasterModelGrid(nrows, ncols, dx)

#create the fields in the grid
mg.create_node_array_zeros('planet_surface__elevation')
z = mg.create_node_array_zeros() + init_elev
mg['node'][ 'planet_surface__elevation'] = z + numpy.random.rand(len(z))/1000.

print( 'Running ...' )

#instantiate the components:
fr = FlowRouter(mg)
sp = StreamPowerEroder(mg, './drive_sp_params.txt')

#perform the loop:
elapsed_time = 0. #total time in simulation
while elapsed_time < time_to_run:
    print elapsed_time
    if elapsed_time+dt<time_to_run:
        dt = time_to_run - elapsed_time
    mg = fr.route_flow(grid=mg)
    mg,_,_ = sp.erode(mg, dt, node_drainage_areas='drainage_area', slopes_at_nodes='steepest_slope')
    
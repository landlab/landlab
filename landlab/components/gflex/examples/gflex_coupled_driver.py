# -*- coding: utf-8 -*-
"""
A driver for our version of AW's gFlex component.

Created on Fri Feb 20 11:17:52 2015

@author: danhobley
"""
from __future__ import print_function

import numpy as np
import pylab

from landlab.components.gFlex.flexure import gFlex
from landlab.components.flow_routing import FlowRouter
from landlab.components.stream_power import StreamPowerEroder, FastscapeEroder
from landlab import RasterModelGrid
from landlab import ModelParameterDictionary
from landlab.plot.imshow import imshow_node_grid

inputs = ModelParameterDictionary('./coupled_SP_gflex_params.txt')
nrows = inputs.read_int('nrows')
ncols = inputs.read_int('ncols')
dx = inputs.read_float('dx')
dt = inputs.read_float('dt')
time_to_run = inputs.read_float('run_time')
init_elev = inputs.read_float('init_elev')
uplift_perstep = inputs.read_float('uplift_rate')*dt
rock_stress_param = inputs.read_float('rock_density')*9.81

mg = RasterModelGrid(nrows, ncols, dx)

#create the fields in the grid
mg.add_zeros('topographic__elevation', at='node')
z = mg.zeros(at='node') + init_elev
mg['node'][ 'topographic__elevation'] = z + np.random.rand(len(z))/1000.

#make some surface load stresses in a field to test
mg.at_node['surface_load__stress'] = np.zeros(nrows*ncols, dtype=float)

#instantiate:
gf = gFlex(mg, './coupled_SP_gflex_params.txt')
fsp = FastscapeEroder(mg, './coupled_SP_gflex_params.txt')
sp = StreamPowerEroder(mg, './coupled_SP_gflex_params.txt')
fr = FlowRouter(mg)

#perform the loop:
elapsed_time = 0. #total time in simulation
while elapsed_time < time_to_run:
    print(elapsed_time)
    if elapsed_time+dt>time_to_run:
        print("Short step!")
        dt = time_to_run - elapsed_time
    mg = fr.route_flow()
    #mg = fsp.erode(mg)
    mg,_,_ = sp.erode(mg, dt, node_drainage_areas='drainage_area', slopes_at_nodes='topographic__steepest_slope')
    mg.at_node['surface_load__stress'] = (mg.at_node['topographic__elevation']+1000)*rock_stress_param
    gf.flex_lithosphere()
    mg.at_node['topographic__elevation'][mg.number_of_nodes//4:3.*mg.number_of_nodes//4] += uplift_perstep
    elapsed_time += dt

pylab.figure(1)
im = imshow_node_grid(mg, 'topographic__elevation')  # display a colored image

pylab.figure(2)
im = imshow_node_grid(mg, 'lithosphere_surface__elevation_increment')

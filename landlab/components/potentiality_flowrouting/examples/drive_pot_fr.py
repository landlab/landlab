# -*- coding: utf-8 -*-
"""drive_pot_fr

Drive the landlab potentiality flow routing component.

Created on Wed Mar 4 2015

@author: danhobley
"""

from landlab import RasterModelGrid, ModelParameterDictionary
from landlab.plot.imshow import imshow_node_grid
import numpy as np
from pylab import imshow, show, contour, figure, clabel, quiver
from landlab.components.potentiality_flowrouting.route_flow_by_boundary import PotentialityFlowRouter
from landlab.components.flow_routing.route_flow_dn import FlowRouter
from landlab.components.stream_power.stream_power import StreamPowerEroder

nrows = 50
ncols = 50

mg = RasterModelGrid(nrows, ncols)
z = (mg.node_y)*0.5
#z = np.sqrt(mg.node_x**2 + mg.node_y**2)
#z = mg.node_x + mg.node_y
#val_to_replace_with = z.reshape((nrows,ncols))[3,3]
#z.reshape((nrows,ncols))[2:5,:] = val_to_replace_with
mg.at_node['topographic_elevation'] = z

mg.set_fixed_value_boundaries_at_grid_edges(True, True, True, True)

mg.at_node['water__volume_flux_in'] = np.ones_like(z)

pfr = PotentialityFlowRouter(mg, 'pot_fr_params.txt')

pfr.route_flow(route_on_diagonals=False)

figure(1)
imshow_node_grid(mg, 'water__volume_flux_magnitude')
figure(2)
imshow_node_grid(mg, 'topographic_elevation')
show()

print mg.at_node['water__volume_flux_magnitude'].reshape((nrows,ncols))

#now a run with a grid...
#make a topo to test on:
mg.at_node['topographic_elevation'] = np.zeros(mg.number_of_nodes)
fr = FlowRouter(mg)
sp = StreamPowerEroder(mg, './pot_fr_params.txt')
inputs = ModelParameterDictionary('./drive_sp_params.txt')
dt = inputs.read_float('dt')
time_to_run = inputs.read_float('run_time')
uplift = inputs.read_float('uplift_rate')
mg.set_fixed_value_boundaries_at_grid_edges(True, True, True, True)
#make some K values in a field to test 
mg.at_node['K_values'] = 0.1+np.random.rand(nrows*ncols)/10.
elapsed_time = 0. #total time in simulation
while elapsed_time < time_to_run:
    #add uplift
    mg.at_node['topographic_elevation'][mg.core_nodes] += uplift*dt
    print elapsed_time
    if elapsed_time+dt>time_to_run:
        print "Short step!"
        dt = time_to_run - elapsed_time
    mg = fr.route_flow(grid=mg)
    mg,_,_ = sp.erode(mg, dt, node_drainage_areas='drainage_area', slopes_at_nodes='steepest_slope', K_if_used='K_values')
    elapsed_time += dt

pfr = PotentialityFlowRouter(mg, 'pot_fr_params.txt')
pfr.route_flow()#route_on_diagonals=False)

figure(1)
imshow_node_grid(mg, 'water__volume_flux_magnitude')
figure(2)
imshow_node_grid(mg, 'drainage_area')
show()
imshow_node_grid(mg, 'potentiality_field')
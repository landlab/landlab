# -*- coding: utf-8 -*-
"""drive_pot_fr

Drive the landlab potentiality flow routing component.

Created on Wed Mar 4 2015

@author: danhobley
"""
from __future__ import print_function

from landlab import RasterModelGrid, ModelParameterDictionary
from landlab.plot.imshow import imshow_node_grid
import numpy as np
from pylab import imshow, show, contour, figure, clabel, quiver, plot, close
from landlab.components.potentiality_flowrouting import PotentialityFlowRouter
from landlab.components.flow_routing import FlowRouter
from landlab.components.stream_power import FastscapeEroder

nrows = 100
ncols = 100
dx = 1000.

close('all')
mg = RasterModelGrid(nrows, ncols, dx)
z = (3000.-mg.node_x)*0.5
#z = -mg.node_x-mg.node_y
#z = np.sqrt(mg.node_x**2 + mg.node_y**2)
#z = mg.node_x + mg.node_y
#val_to_replace_with = z.reshape((nrows,ncols))[3,3]
#z.reshape((nrows,ncols))[2:5,:] = val_to_replace_with
mg.at_node['topographic__elevation'] = z

#mg.set_fixed_value_boundaries_at_grid_edges(True, True, True, True)
mg.set_closed_boundaries_at_grid_edges(False, True, True, True)
figure(3)
imshow_node_grid(mg, mg.status_at_node)

mg.at_node['water__unit_flux_in'] = np.ones_like(z)

pfr = PotentialityFlowRouter(mg, 'pot_fr_params.txt')

pfr.route_flow(route_on_diagonals=True)

figure(1)
imshow_node_grid(mg, 'surface_water__discharge')
figure(2)
imshow_node_grid(mg, 'topographic__elevation')

out_sum = np.sum(mg.at_node['surface_water__discharge'].reshape((nrows,ncols))[-3,:])
print(out_sum)
print(np.sum(mg.at_node['water__unit_flux_in']), np.sum(mg.at_node['surface_water__discharge'][mg.boundary_nodes]))
print(out_sum - np.sum(mg.at_node['water__unit_flux_in']))

show()

print(mg.at_node['surface_water__discharge'].reshape((nrows,ncols)))
print(np.sum(mg.at_node['surface_water__discharge'].reshape((nrows,ncols)),axis=0)/3.)

#now a run with a grid...

#mg = RasterModelGrid(nrows, ncols, dx)
##make a topo to test on:
#mg.at_node['topographic__elevation'] = np.zeros(mg.number_of_nodes)
#fr = FlowRouter(mg)
#fsp = FastscapeEroder(mg, './pot_fr_params.txt')
#inputs = ModelParameterDictionary('./pot_fr_params.txt')
#dt = inputs.read_float('dt')
#time_to_run = inputs.read_float('run_time')
#uplift = inputs.read_float('uplift_rate')
#mg.set_fixed_value_boundaries_at_grid_edges(True, True, True, True)
##make some K values in a field to test
#mg.at_node['K_values'] = 0.1+np.random.rand(nrows*ncols)/10.
##elapsed_time = 0. #total time in simulation
##while elapsed_time < time_to_run:
##    #add uplift
##    mg.at_node['topographic__elevation'][mg.core_nodes] += uplift*dt
##    print elapsed_time
##    if elapsed_time+dt>time_to_run:
##        print "Short step!"
##        dt = time_to_run - elapsed_time
##    mg = fr.route_flow(grid=mg)
##    mg,_,_ = sp.erode(mg)
##    elapsed_time += dt

inputs = ModelParameterDictionary('./pot_fr_params.txt')
nrows = 200#inputs.read_int('nrows')
ncols = 200#inputs.read_int('ncols')
dx = inputs.read_float('dx')
dt = inputs.read_float('dt')
time_to_run = inputs.read_float('run_time')
#nt needs defining
uplift = inputs.read_float('uplift_rate')
init_elev = inputs.read_float('init_elev')

mg = RasterModelGrid(nrows, ncols, dx)

#create the fields in the grid
mg.add_zeros('topographic__elevation', at='node')
z = mg.zeros(at='node') + init_elev
mg['node']['topographic__elevation'] = z + np.random.rand(len(z))/1000.

#make some K values in a field to test
mg.at_node['K_values'] = 0.00001+np.random.rand(nrows*ncols)/100000.

#mg.at_node['water__unit_flux_in'] = dx*dx*np.ones_like(z)
mg.at_node['water__unit_flux_in'] = dx*dx*np.ones_like(z)*100./(60.*60.*24.*365.25) #remember, flux is /sec, so this is a small number!
#mg.set_closed_boundaries_at_grid_edges(False, False, True, True)
#mg.set_closed_boundaries_at_grid_edges(True, False, True, True)

print( 'Running ...' )

#instantiate the components:
fr = FlowRouter(mg)
#load the Fastscape module too, to allow direct comparison
fsp = FastscapeEroder(mg, './pot_fr_params.txt')

#perform the loop:
elapsed_time = 0. #total time in simulation
while elapsed_time < time_to_run:
    print(elapsed_time)
    if elapsed_time+dt>time_to_run:
        print("Short step!")
        dt = time_to_run - elapsed_time
    mg = fr.route_flow()
    #print 'Area: ', numpy.max(mg.at_node['drainage_area'])
    #mg = fsp.erode(mg)
    mg = fsp.erode(mg, K_if_used='K_values')
    #mg,_,_ = sp.erode(mg, dt, node_drainage_areas='drainage_area', slopes_at_nodes='topographic__steepest_slope')
    #add uplift
    mg.at_node['topographic__elevation'][mg.core_nodes] += uplift*dt
    elapsed_time += dt

pfr = PotentialityFlowRouter(mg, 'pot_fr_params.txt')
pfr.route_flow(return_components=True)#route_on_diagonals=False)

figure('Topo')
imshow_node_grid(mg, 'topographic__elevation')
figure('Potentiality flow fluxes')
imshow_node_grid(mg, 'surface_water__discharge')
figure('D8 drainage areas')
imshow_node_grid(mg, 'drainage_area')
figure('K (core only)')
mg.at_node['flow__potential'][mg.boundary_nodes] = 0.
imshow_node_grid(mg, 'flow__potential')
depths = np.where(mg.at_node['surface_water__depth']>1.e6,0.,mg.at_node['surface_water__depth'])
figure('depths')
imshow_node_grid(mg, depths, var_name='depths (Manning)')
figure('xcomp')
imshow_node_grid(mg, 'water__discharge_x_component')
figure('ycomp')
imshow_node_grid(mg, 'water__discharge_y_component')

print('flux in per node: ', mg.at_node['water__unit_flux_in'][0])
print('water in, water out: ', np.sum(mg.at_node['water__unit_flux_in'][mg.core_nodes]), np.sum(mg.at_node['surface_water__discharge'][mg.boundary_nodes]))
print(-np.sum(mg.at_node['surface_water__discharge'][mg.boundary_nodes]) + np.sum(mg.at_node['water__unit_flux_in'][mg.core_nodes]))
print((-np.sum(mg.at_node['surface_water__discharge'][mg.boundary_nodes]) + np.sum(mg.at_node['water__unit_flux_in'][mg.core_nodes]))/mg.at_node['water__unit_flux_in'][0])

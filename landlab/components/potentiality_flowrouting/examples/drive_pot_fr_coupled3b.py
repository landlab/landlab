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
from landlab.components.potentiality_flowrouting.route_flow_by_boundary import PotentialityFlowRouter
from landlab.components.flow_routing.route_flow_dn import FlowRouter
from landlab.components.stream_power.fastscape_stream_power import SPEroder
from landlab.grid.mappers import map_link_end_node_max_value_to_link

inputs = ModelParameterDictionary('./pot_fr_params.txt')
nrows = 50#inputs.read_int('nrows')
ncols = 50#inputs.read_int('ncols')
dx = inputs.read_float('dx')
init_elev = inputs.read_float('init_elev')

mg = RasterModelGrid(nrows, ncols, dx)

# attempt to implement diffusion with flow routing...

#modify the fields in the grid
z = mg.create_node_array_zeros() + init_elev
z_slope = (49000. - mg.node_y)/mg.node_y.max()/20.
mg.at_node['topographic__elevation'] = z + z_slope #+ np.random.rand(len(z))/1000.
mg.create_node_array_zeros('water__volume_flux_in')

#Set boundary conditions
inlet_node = np.array((int((1.5*mg.number_of_node_columns)//1)))
section_col = int((0.5*mg.number_of_node_columns)//1)
mg.at_node['topographic__elevation'][section_col] = 1.
mg.set_closed_boundaries_at_grid_edges(True, False, False, False)
mg.set_fixed_value_boundaries_at_grid_edges(False, True, True, True)
mg.node_status[section_col] = 2
mg.update_links_nodes_cells_to_new_BCs()
mg.at_node['water__volume_flux_in'].fill(0.)
mg.at_node['water__volume_flux_in'][inlet_node] = 1.
pfr = PotentialityFlowRouter(mg, 'pot_fr_params.txt')

interior_nodes = mg.core_nodes

#store profiles here
section_downfan = []

# do the loop
for i in xrange(3000):
    #mg.at_node['topographic__elevation'][inlet_node] = 1.
    #maintain flux like this now instead:
    mg.at_node['topographic__elevation'][section_col] = mg.at_node['topographic__elevation'][inlet_node]+1.
    pfr.route_flow(route_on_diagonals=True)
    #imshow(mg, 'water__volume_flux_magnitude')
    #show()
    kd = mg.at_node['water__volume_flux_magnitude']   # 0.01 m2 per year
    # dt = np.nanmin(0.2*mg.dx*mg.dx/kd)   # CFL condition
    dt = 0.5
    g = mg.calculate_gradients_at_active_links(mg.at_node['topographic__elevation'])
    map_link_end_node_max_value_to_link(mg, 'water__volume_flux_magnitude')
    kd_link = 1.e6*mg.at_link['water__volume_flux_magnitude'][mg.active_links]
    qs = -kd_link*g
    dqsdx = mg.calculate_flux_divergence_at_nodes(qs)
    dzdt = -dqsdx
    mg.at_node['topographic__elevation'][interior_nodes] += dzdt[interior_nodes]*dt
    if i%50==0:
        print('loop '+str(i))
        section_downfan.append(mg.node_vector_to_raster(mg.at_node['topographic__elevation'])[1:,section_col].copy())

#drop the BL HARD
mg.at_node['topographic__elevation'][mg.top_edge_node_ids()[mg.number_of_node_columns//2]] =- 50.

for i in xrange(3000):
    #mg.at_node['topographic__elevation'][inlet_node] = 1.
    #maintain flux like this now instead:
    mg.at_node['topographic__elevation'][section_col] = mg.at_node['topographic__elevation'][inlet_node]+1.
    pfr.route_flow(route_on_diagonals=True)
    #imshow(mg, 'water__volume_flux_magnitude')
    #show()
    kd = mg.at_node['water__volume_flux_magnitude']   # 0.01 m2 per year
    # dt = np.nanmin(0.2*mg.dx*mg.dx/kd)   # CFL condition
    dt = 0.5
    g = mg.calculate_gradients_at_active_links(mg.at_node['topographic__elevation'])
    map_link_end_node_max_value_to_link(mg, 'water__volume_flux_magnitude')
    kd_link = 1.e6*mg.at_link['water__volume_flux_magnitude'][mg.active_links]
    qs = -kd_link*g
    dqsdx = mg.calculate_flux_divergence_at_nodes(qs)
    dzdt = -dqsdx
    mg.at_node['topographic__elevation'][interior_nodes] += dzdt[interior_nodes]*dt
    if i%50==0:
        print('loop '+str(i))
        section_downfan.append(mg.node_vector_to_raster(mg.at_node['topographic__elevation'])[1:,section_col].copy())

figure(1)
imshow_node_grid(mg, 'topographic__elevation')
figure(2)
imshow_node_grid(mg, mg.hillshade(), cmap='bone')
figure(3)
imshow_node_grid(mg, 'water__volume_flux_magnitude', cmap='Blues_r')
figure(4)
for i in xrange(len(section_downfan)):
    plot(section_downfan[i], '-')

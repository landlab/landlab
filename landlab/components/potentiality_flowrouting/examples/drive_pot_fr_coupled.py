# -*- coding: utf-8 -*-
"""drive_pot_fr

Drive the landlab potentiality flow routing component.

Created on Wed Mar 4 2015

@author: danhobley
"""
from __future__ import print_function

from six.moves import range

from landlab import RasterModelGrid, ModelParameterDictionary
from landlab.plot.imshow import imshow_node_grid
import numpy as np
from pylab import imshow, show, contour, figure, clabel, quiver, plot, close
from landlab.components.potentiality_flowrouting import PotentialityFlowRouter
from landlab.components.flow_routing import FlowRouter
from landlab.components.stream_power import FastscapeEroder
from landlab.grid.mappers import map_link_end_node_max_value_to_link

inputs = ModelParameterDictionary('./pot_fr_params.txt')
nrows = 50#inputs.read_int('nrows')
ncols = 50#inputs.read_int('ncols')
dx = inputs.read_float('dx')
init_elev = inputs.read_float('init_elev')

mg = RasterModelGrid(nrows, ncols, dx)

# attempt to implement diffusion with flow routing...

#modify the fields in the grid
z = mg.zeros(at='node') + init_elev
mg.at_node['topographic__elevation'] = z + np.random.rand(len(z))/1000.
mg.add_zeros('water__unit_flux_in', at='node')

#Set boundary conditions
mg.set_closed_boundaries_at_grid_edges(False, True, True, True)
mg.set_fixed_value_boundaries_at_grid_edges(False, False, False, True)
inlet_node = np.array((mg.number_of_node_columns + 1))
mg.at_node['water__unit_flux_in'].fill(0.)
mg.at_node['water__unit_flux_in'][inlet_node] = 1.
pfr = PotentialityFlowRouter(mg, 'pot_fr_params.txt')

interior_nodes = mg.core_nodes

# do the loop
for i in range(2000):
    if i%50==0:
        print('loop '+str(i))
    mg.at_node['topographic__elevation'][inlet_node] = 1.
    pfr.route_flow(route_on_diagonals=True)
    #imshow(mg, 'surface_water__discharge')
    #show()
    kd = mg.at_node['surface_water__discharge']   # 0.01 m2 per year
    # dt = np.nanmin(0.2*mg.dx*mg.dx/kd)   # CFL condition
    dt = 0.5
    g = mg.calc_grad_of_active_link(mg.at_node['topographic__elevation'])
    map_link_end_node_max_value_to_link(mg, 'surface_water__discharge')
    kd_link = 1.e6*mg.at_link['surface_water__discharge'][mg.active_links]
    qs = -kd_link*g
    dqsdx = mg.calculate_flux_divergence_at_nodes(qs)
    dzdt = -dqsdx
    mg.at_node['topographic__elevation'][interior_nodes] += dzdt[interior_nodes]*dt

figure(1)
imshow_node_grid(mg, 'topographic__elevation')
figure(2)
imshow_node_grid(mg, 'surface_water__depth')
figure(3)
imshow_node_grid(mg, 'surface_water__discharge')

# -*- coding: utf-8 -*-
"""
Created on Wed May 22 13:50:43 2019

@author: abby
"""

"""
Created on Thu Apr  4 11:57:07 2019

@author: abby

First try of a lateral erosion test after many years.....

Modification of K1A2spin30_lat.py
"""

from IPython import get_ipython
ipython = get_ipython()
ipython.magic("reset -f")

import os
import numpy as np
from pylab import *
from landlab.plot.imshow import imshow_grid
from landlab.components import FlowAccumulator

from landlab import RasterModelGrid, CLOSED_BOUNDARY, FIXED_VALUE_BOUNDARY
from matplotlib.pyplot import figure, show, plot, xlabel, ylabel, title

from numpy import testing

from random import uniform
from landlab.components import LateralEroder
from landlab import load_params

from landlab.utils import structured_grid

import time


#%% 
"""
test to see if sediment flux of model matches SS analytical solution
"""
tic=time.time()

#some helper and parameter values
U = 0.0005 #m/yr
dt = 100 #years

nr = 5
nc = 5
dx=10
#instantiate grid
mg = RasterModelGrid(nr, nc, dx)

for edge in (mg.nodes_at_top_edge, mg.nodes_at_bottom_edge, mg.nodes_at_left_edge, mg.nodes_at_right_edge):
    mg.status_at_node[edge] = CLOSED_BOUNDARY
for edge in (mg.nodes_at_bottom_edge):
    mg.status_at_node[edge] = FIXED_VALUE_BOUNDARY

z = mg.add_zeros('node', 'topographic__elevation')
ir2 = np.random.uniform(low=0.0, high=0.5, size=(z.size))
loading_vector = np.linspace(1,2.5,num=nr)
ramp = np.repeat(loading_vector, nc)
ramp += ir2

#z += initial_roughness
z += ramp
fa = FlowAccumulator(mg,
                     surface='topographic__elevation',
                     flow_director='FlowDirectorD8',
                     runoff_rate=None,
                     depression_finder=None)#"DepressionFinderAndRouter", router="D8")

latero = LateralEroder(mg,latero_mech="UC", Kv=0.0001, Kl_ratio=1.5, solver="basic")

for i in range(1000):
    fa.run_one_step()    #flow accumulator
#    (da, q) = fa.accumulate_flow()
    (mg, dzlat)=latero.run_one_step(mg,dt,)
    mg.at_node['topographic__elevation'][mg.core_nodes] += U*dt    #uplift the landscape


    # compare numerical and analytical sediment flux solutions
    
    ####****NOTE. MAY 22 2019: THIS WON'T WORK. I need to have a qs variable as well as qsin.
    # have now a qs in model. it works.
num_sedflux = mg.at_node["qs"][mg.core_nodes]
analytical_sedflux = U * mg.at_node["drainage_area"][mg.core_nodes] 
print(num_sedflux)
print(analytical_sedflux)

toc=time.time()
print ("total time", toc-tic)

from landlab.plot.drainage_plot import drainage_plot
plt.figure()
drainage_plot(mg)
# test for match with anakytical sediment flux. note tha they are off a little 
#because of the lateral erosion
testing.assert_array_almost_equal(
    num_sedflux,
    analytical_sedflux,
    decimal=2,
    err_msg="SPACE transport-limited sediment flux test failed",
    verbose=True,
)


#%%
"""
variable bedrock erodibility
"""

U = 0.005 #m/yr
dt = 10 #years

nr = 5
nc = 5
nnodes = nr*nc
dx=10
#instantiate grid
mg = RasterModelGrid(nr, nc, dx)
#below gives you node id of boundary nodes
boundary_nodes=structured_grid.boundary_nodes((nr,nc))
#rg.set_closed_boundaries_at_grid_edges(True, True, True, False)    #bottom is open
for edge in (mg.nodes_at_top_edge, mg.nodes_at_bottom_edge, mg.nodes_at_left_edge, mg.nodes_at_right_edge):
    mg.status_at_node[edge] = CLOSED_BOUNDARY
for edge in (mg.nodes_at_bottom_edge):
    mg.status_at_node[edge] = FIXED_VALUE_BOUNDARY

z = mg.add_zeros('node', 'topographic__elevation')
# add some roughness, as this lets "natural" channel planforms arise
    
loading_vector = np.linspace(1,4,num=nr)
ramp = np.repeat(loading_vector, nc)
ramp += np.random.random_sample(nnodes)*0.8

z += ramp

Kvar=0.006*np.ones(nr*nc)
Kvar[0:9]=0.06
#set up fluvial lateral erosion incision component
fa = FlowAccumulator(mg,
                     surface='topographic__elevation',
                     flow_director='FlowDirectorD8',
                     runoff_rate=None,
                     depression_finder=None)#"DepressionFinderAndRouter", router="D8")

fa.run_one_step()
#***NOTE, YOU MUST USE ADAPTIVE TIME STEPPER FOR variable K, or you may get strange 
# topography
latero = LateralEroder(mg,latero_mech="UC", Kv=Kvar, solver="adaptive")

for i in range(50):

    fa.run_one_step()    #flow accumulator
    #erode the landscape with lateral erosion
    (mg, dzlat)=latero.run_one_step(mg,dt,)
    mg.at_node['topographic__elevation'][mg.core_nodes] += U*dt 


#***NOTE: have to make a test with some analytical solution with K and slope.
print(mg.at_node['topographic__steepest_slope'][mg.core_nodes])    
    
from landlab.plot.drainage_plot import drainage_plot
plt.figure()
drainage_plot(mg)

figure()
im = imshow_grid(mg, Kvar, grid_units = ['m','m'],
             var_name='kv', )

#%%
"""
steady inlet

NOTE: HOW do I test if this is working properly? just ask if da[node]=set value?
NOTE:******* CHECK if the sediment and da are matched.
"""


U = 0.005 #m/yr
dt = 10 #years

nr = 5
nc = 5
nnodes = nr*nc
dx=10
#instantiate grid
mg = RasterModelGrid(nr, nc, dx)
#below gives you node id of boundary nodes
boundary_nodes=structured_grid.boundary_nodes((nr,nc))
#rg.set_closed_boundaries_at_grid_edges(True, True, True, False)    #bottom is open
for edge in (mg.nodes_at_top_edge, mg.nodes_at_bottom_edge, mg.nodes_at_left_edge, mg.nodes_at_right_edge):
    mg.status_at_node[edge] = CLOSED_BOUNDARY
for edge in (mg.nodes_at_bottom_edge):
    mg.status_at_node[edge] = FIXED_VALUE_BOUNDARY

z = mg.add_zeros('node', 'topographic__elevation')
# add some roughness, as this lets "natural" channel planforms arise
    
loading_vector = np.linspace(1,2.5,num=nr)
ramp = np.repeat(loading_vector, nc)
ramp += np.random.random_sample(nnodes)*0.8
z += ramp


latero = LateralEroder(mg,latero_mech="UC", Kv=0.0001, Kl_ratio=1.5, inlet_on=True, 
                       inlet_node=17, inlet_area=1e3, qsinlet=5.)

fa = FlowAccumulator(mg,
                     surface='topographic__elevation',
                     flow_director='FlowDirectorD8',
                     runoff_rate=None,
                     depression_finder=None)#"DepressionFinderAndRouter", router="D8")
for i in range(100):

#    fa.run_one_step()    #flow accumulator
    fa.accumulate_flow()
    #erode the landscape with lateral erosion
    # if you want the lateral erosion component to see the qsin you've added, must provide inlet node and qsinlet
    (mg, dzlat)=latero.run_one_step(mg,dt,)
    mg.at_node['topographic__elevation'][mg.core_nodes] += U*dt    #uplift the landscape
    
plt.figure()
drainage_plot(mg)
da=mg.at_node['surface_water__discharge']/dx**2
figure()
im = imshow_grid(mg, 'surface_water__discharge', cmap="jet",grid_units = ['m','m'],
             var_name='surface water discharge (m)')
figure()
im = imshow_grid(mg, 'qs', cmap="bone",vmin=2, vmax=6,grid_units = ['m','m'])
figure()
im = imshow_grid(mg, da, cmap="jet",grid_units = ['m','m'])
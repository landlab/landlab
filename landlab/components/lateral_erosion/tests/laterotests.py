# -*- coding: utf-8 -*-
"""
Created on Wed May 22 13:50:43 2019

@author: abby
"""

"""
Created on Thu Apr  4 11:57:07 2019

@author: abby

Doc tests and unit tests for lateral erosion.
"""

from IPython import get_ipython
ipython = get_ipython()
ipython.magic("reset -f")

import numpy as np
from landlab import RasterModelGrid, CLOSED_BOUNDARY, FIXED_VALUE_BOUNDARY
from matplotlib.pyplot import figure, show, plot, xlabel, ylabel, title

from numpy import testing
from random import uniform
from landlab import load_params

from landlab.utils import structured_grid

import time

from landlab.plot.imshow import imshow_grid
from matplotlib.pyplot import figure, show, plot, xlabel, ylabel, title
from landlab.components import FlowAccumulator, LateralEroder
#%% test that sets up a simple, pre-defined drainage network and compares
# the lateral node that is eroded, the volume of lateral eorsion, and the elevation
# of the landscape after one timestep
nr = 5
nc = 5
nnodes = nr*nc
dx=1
#instantiate grid
mg = RasterModelGrid(nr, nc, dx)

#rg.set_closed_boundaries_at_grid_edges(True, True, True, False)    #bottom is open
for edge in (mg.nodes_at_top_edge, mg.nodes_at_bottom_edge, mg.nodes_at_left_edge, mg.nodes_at_right_edge):
    mg.status_at_node[edge] = CLOSED_BOUNDARY
for edge in (mg.nodes_at_bottom_edge):
    mg.status_at_node[edge] = FIXED_VALUE_BOUNDARY

z = mg.add_zeros('node', 'topographic__elevation')
loading_vector = np.linspace(1,4,num=nr)
ramp = np.repeat(loading_vector, nc)
#the tweaks to elevation below make lateral node at node 7
z += ramp
z[11]-=0.9
z[12]-=0.4
z[8]-=0.001
#print("zbefore", mg.at_node['topographic__elevation'].reshape(5,5))
fa = FlowAccumulator(mg,
                     surface='topographic__elevation',
                     flow_director='FlowDirectorD8',
                     runoff_rate=None,
                     depression_finder=None)#"DepressionFinderAndRouter", router="D8")
latero = LateralEroder(mg,latero_mech="UC", Kv=0.1, Kl_ratio=1.5)
fa.accumulate_flow()
    #erode the landscape with lateral erosion
(mg, dzlat,)=latero.run_one_step(mg,dt=1.,)

#print("zafter", mg.at_node['topographic__elevation'].reshape(5,5))
qsname=mg['node'][ 'qs_in']
#print("qsin", qsname.reshape(5,5))
vlname=mg['node'][ 'volume__lateral_erosion']
#print("vollatout", vlname)
#print('vollat[7]', vlname[7])
pred_vollat=0.00045158164
pred_zafter=np.array([1., 1., 1., 1., 1., 1.75, 1.675, 1.675, 1.6731154, 1.75, 2.5,  1.66418779,  2.06181623,  2.4249,  2.5, 3.25,  3.085,  3.13332738,  3.16868272,  3.25, 4., 4., 4., 4., 4.])
#ztf=np.allclose(mg.at_node['topographic__elevation'], pred_zafter)
#tf=np.allclose(vlname[7], pred_vollat)

testing.assert_array_almost_equal(
    mg.at_node['topographic__elevation'],
    pred_zafter,
    decimal=8,
    err_msg="LatEro basic erosion test failed",
    verbose=True,
)
testing.assert_array_almost_equal(
    vlname[7],
    pred_vollat,
    decimal=8,
    err_msg="LatEro volume lateral erosion failed",
    verbose=True,
)
#from landlab.plot.drainage_plot import drainage_plot
#figure()
#drainage_plot(mg)
#print(delta)
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
num_sedflux = mg.at_node["qs"][mg.core_nodes]
analytical_sedflux = U * mg.at_node["drainage_area"][mg.core_nodes] 

# test for match with anakytical sediment flux. note tha they are off a little 
#because of the lateral erosion
testing.assert_array_almost_equal(
    num_sedflux,
    analytical_sedflux,
    decimal=2,
    err_msg="LatEro transport-limited sediment flux test failed",
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
for edge in (mg.nodes_at_top_edge, mg.nodes_at_bottom_edge, mg.nodes_at_left_edge, mg.nodes_at_right_edge):
    mg.status_at_node[edge] = CLOSED_BOUNDARY
for edge in (mg.nodes_at_bottom_edge):
    mg.status_at_node[edge] = FIXED_VALUE_BOUNDARY

z = mg.add_zeros('node', 'topographic__elevation')
loading_vector = np.linspace(1,4,num=nr)
ramp = np.repeat(loading_vector, nc)
ramp += np.random.random_sample(nnodes)*0.8
z += ramp

#below, array of different bedrock erodibilities
Kvar=0.006*np.ones(nr*nc)
Kvar[0:9]=0.06
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
    
#from landlab.plot.drainage_plot import drainage_plot
#plt.figure()
#drainage_plot(mg)
#
#figure()
#im = imshow_grid(mg, Kvar, grid_units = ['m','m'],
#             var_name='kv', )

#%%
"""
steady inlet

"""

U = 0.005 #m/yr
dt = 10 #years

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
loading_vector = np.linspace(1,2.5,num=nr)
ramp = np.repeat(loading_vector, nc)
ramp += np.random.random_sample(mg.number_of_nodes)*0.8
z += ramp

latero = LateralEroder(mg,latero_mech="UC", Kv=0.001, Kl_ratio=1.5, inlet_on=True, 
                       inlet_node=17, inlet_area=500, qsinlet=2.5)

fa = FlowAccumulator(mg,
                     surface='topographic__elevation',
                     flow_director='FlowDirectorD8',
                     runoff_rate=None,
                     depression_finder=None)#"DepressionFinderAndRouter", router="D8")
for i in range(1000):

    fa.run_one_step()    #flow accumulator
    #erode the landscape with lateral erosion
    (mg, dzlat)=latero.run_one_step(mg,dt,)
    mg.at_node['topographic__elevation'][mg.core_nodes] += U*dt    #uplift the landscape

da=mg.at_node['surface_water__discharge']/dx**2
num_sedflux = mg.at_node["qs"]
analytical_sedflux = U * da 

# test for match with analytical sediment flux. note that the values are off a little 
#because of the lateral erosion
testing.assert_array_almost_equal(
    num_sedflux,
    analytical_sedflux,
    decimal=2,
    err_msg="LatEro inlet transport-limited sediment flux test failed",
    verbose=True,
)
    
#%%
"""
time varying inlet

"""

U = 0.005 #m/yr
dt = 10 #years

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
loading_vector = np.linspace(1,2.5,num=nr)
ramp = np.repeat(loading_vector, nc)
ramp += np.random.random_sample(mg.number_of_nodes)*0.8
z += ramp

latero = LateralEroder(mg,latero_mech="UC", Kv=0.001, Kl_ratio=1.5, inlet_on=True, 
                       inlet_node=17, inlet_area=500, qsinlet=2.5)

fa = FlowAccumulator(mg,
                     surface='topographic__elevation',
                     flow_director='FlowDirectorD8',
                     runoff_rate=None,
                     depression_finder=None)#"DepressionFinderAndRouter", router="D8")
inareats=None
sedints=None
for i in range(1000):
    if i == 200:
        #when i equals 200, increase drainage area and sediment flux at inlet
        inareats = 1000
        sedints = 5
    fa.run_one_step()    #flow accumulator
    (mg, dzlat)=latero.run_one_step(mg,dt,inlet_area_ts=inareats, qsinlet_ts=sedints)
    mg.at_node['topographic__elevation'][mg.core_nodes] += U*dt    #uplift the landscape

da=mg.at_node['surface_water__discharge']/dx**2
num_sedflux = mg.at_node["qs"]
analytical_sedflux = U * da 

# test for match with analytical sediment flux. note that the values are off a little 
#because of the lateral erosion
testing.assert_array_almost_equal(
    num_sedflux,
    analytical_sedflux,
    decimal=1,
    err_msg="LatEro time varying inlet transport-limited sediment flux test failed",
    verbose=True,
)
toc=time.time()
print ("total time", toc-tic)
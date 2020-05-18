#!/usr/bin/env python
# coding: utf-8

# # Landlab Synthetic NetCDF Builder
# 
# <img src="https://www.washington.edu/brand/files/2014/09/W-Logo_Purple_Hex.png" style="float:right;width:200px;padding:20px">   
# 
# 
# <br />
# This Jupyter Notebook runs the Landlab LandslideProbability component on a synthetic 
# Landlab grid using four depth to water table options to replace recharge options described in the paper: <br />
# #### Strauch et al. 2018. A hydro-climatological approach to predicting regional landslide probability using Landlab. Earth Surface Dynamics, 6, 1-26. <br /> 
# This notebook performs the following functions:<br >
# * Import libraries and set HydroShare variables<br />
# * Create a grid and data fields used to calculate landslide probability<br />
# * Specify Depth to Water Table Distributions to compare four options<br /> 
# * Run LandslideProbability function from Landlab landslide component<br /> 
# * Compare the sensitivity based on four Depth to Water Table options<br /> 
# 

# In[5]:


from landlab import RasterModelGrid
import numpy as np
from landlab.io.netcdf import read_netcdf
from landlab.io.netcdf import write_netcdf


# In[6]:


grid = RasterModelGrid((5, 4))


# In[3]:


grid.nodes


# In[4]:


grid.shape == (5, 4)


# In[23]:


n=50
gridnodes = grid.number_of_nodes
grid_size = grid.number_of_nodes

Demin_value = 2 
Demax_value = 5
distribution1 = 'uniform'
depth_dist = np.random.uniform(Demin_value, Demax_value,size=n)
print('Depth to water table distribution')
print(depth_dist)

mean_depth=np.mean(depth_dist)
grid['node']['soil__mean_watertable_depth']=mean_depth* np.ones(gridnodes)

print('Mean depth to water table from uniform distribution')
print(mean_depth)
print('Mean Depth to water table - uniform for all nodes')
print(grid['node']['soil__mean_watertable_depth'])


# In[24]:


grid.dy, grid.dx


# In[25]:


list(grid.at_node.keys())


# In[26]:


grid.at_node['soil__mean_watertable_depth']


# In[27]:


write_netcdf('synthetic_depth.nc', grid, format='NETCDF3_64BIT', names='soil__mean_watertable_depth')


# In[ ]:





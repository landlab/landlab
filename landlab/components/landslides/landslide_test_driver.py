# -*- coding: utf-8 -*-
"""
Testing driver for landslide.py

Created on Mon Mar 27 10:33:49 2017

@author: rstrauch
"""

# Import libraries
import numpy as np
import os
import matplotlib.pyplot as plt
import matplotlib as mpl
from landlab import RasterModelGrid
from landlab.components.landslides import LandslideProbability

#Create a grid on which to calculate landslide probability.

grid = RasterModelGrid((5, 4), spacing=(0.2, 0.2))
# Add required fields for component.
grid['node']['slope'] = np.random.rand(5,4)
scatter_dat = np.random.random_integers(1, 10, (5,4))
grid['node']['contributing_area']= \
         np.sort(np.random.random_integers(30, 900, (5,4)))
grid['node']['soil_transmissivity']= \
         np.sort(np.random.random_integers(5, 20, (5,4)),-1)
grid['node']['combined_cohesion__mode']= \
         np.sort(np.random.random_integers(30, 900, (5,4)))
grid['node']['combined_cohesion__minimum']= \
         grid.at_node['combined_cohesion__mode'] - scatter_dat
grid['node']['combined_cohesion__maximum']= \
         grid.at_node['combined_cohesion__mode'] + scatter_dat
grid['node']['soil_internal_angle_friction__mode']= \
         np.sort(np.random.random_integers(26, 40, (5,4)))
grid['node']['soil_thickness__mode']= \
         np.sort(np.random.random_integers(1, 10, (5,4)))
grid['node']['soil_density']= \
         2000. * np.ones(grid.number_of_nodes)
         
grid_size = grid.number_of_core_nodes
# number of iterations to run Monte Carlo simulation
n = 10

# Recharge options
#Option 1 - uniform distribution
distribution = 'uniform'
Remin_value = 20.
Remax_value = 120.
#Opton 2 - Lognormal-Uniform
distribution = 'lognormal_uniform'
Remean = 4.
Restandard_deviation = 0.25
#Option 3 - Lognormal-Spatial
distribution = 'lognormal_spatial'
Remean = np.random.random_integers(2,7,grid_size)
Restandard_deviation = np.random.rand(grid_size)
#Option 4 - Fully distributed
VIC_dict = {}
vic_id_dict = {}
fract_dict = {}
# VIC id and recharge array
for vkey in range(2,8):
    VIC_dict[vkey] = np.random.random_integers(20,120,10)
# node id and vic grid ids
for ckey in grid.core_nodes:
    vic_id_dict[ckey] = np.random.random_integers(2,8,2)
# node id and vic grid fraction
for ckey in grid.core_nodes:
    fract_dict[ckey] =  np.random.rand(2)
distribution = 'VIC'
vic_inputs = [VIC_dict,vic_id_dict, fract_dict]

# Instantiate the 'LandslideProbability' component to work on this grid,
# and run it.
#Uniform
LS_prob1 = LandslideProbability(grid,number_of_simulations=n,
    groudwater__recharge_distribution=distribution,
    groundwater__recharge_min_value=Remin_value,
    groundwater__recharge_max_value=Remax_value)
LS_prob1.calculate_landslide_probability()
#lognormal_uniform
LS_prob2 = LandslideProbability(grid,number_of_simulations=n,
    groudwater__recharge_distribution=distribution,
    groundwater__recharge_mean=Remean,
    groundwater__recharge_standard_deviation=Restandard_deviation)
LS_prob2.calculate_landslide_probability()
#lognormal_spatial
LS_prob3 = LandslideProbability(grid,number_of_simulations=n,
    groudwater__recharge_distribution=distribution,
    groundwater__recharge_mean=Remean,
    groundwater__recharge_standard_deviation=Restandard_deviation)
LS_prob3.calculate_landslide_probability()
#VIC
LS_prob4 = LandslideProbability(grid,number_of_simulations=n,
    groudwater__recharge_vic_inputs=vic_inputs)
LS_prob4.calculate_landslide_probability()



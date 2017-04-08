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
gridnum = grid.number_of_nodes
grid['node']['topographic__slope'] = np.random.rand(gridnum)
scatter_dat = np.random.randint(1, 10, gridnum)
grid['node']['topographic__specific_contributing_area']= \
         np.sort(np.random.randint(30, 900, gridnum))
grid['node']['soil__transmissivity']= \
         np.sort(np.random.randint(5, 20, gridnum),-1)
grid['node']['soil__mode_total_cohesion']= \
         np.sort(np.random.randint(30, 900, gridnum))
grid['node']['soil__minimum_total_cohesion']= \
         grid.at_node['soil__mode_total_cohesion'] - scatter_dat
grid['node']['soil__maximum_total_cohesion']= \
         grid.at_node['soil__mode_total_cohesion'] + scatter_dat
grid['node']['soil__internal_friction_angle']= \
         np.sort(np.random.randint(26, 40, gridnum))
grid['node']['soil__thickness']= \
         np.sort(np.random.randint(1, 10, gridnum))
grid['node']['soil__density']= \
         2000. * np.ones(grid.number_of_nodes)
         
grid_size = grid.number_of_core_nodes
# number of iterations to run Monte Carlo simulation
n = 10

# Recharge options
#Option 1 - uniform distribution
distribution1 = 'uniform'
Remin_value = 20.
Remax_value = 120.
#Option 2 - Lognormal-Uniform
distribution2 = 'lognormal_uniform'
Remean = 4.
Restandard_deviation = 0.25
#Option 3 - Lognormal-Spatial
distribution3 = 'lognormal_spatial'
Remean3 = np.random.randint(20,120,grid_size)
Restandard_deviation3 = np.random.rand(grid_size)
#Option 4 - Fully distributed
VIC_dict = {}
vic_id_dict = {}
fract_dict = {}
# VIC id and recharge array
for vkey in range(2,8):
    VIC_dict[vkey] = np.random.randint(20,120,10)
# node id and vic grid ids
for ckey in grid.core_nodes:
    vic_id_dict[ckey] = np.random.randint(2,8,2)
# node id and vic grid fraction
for ckey in grid.core_nodes:
    fract_dict[ckey] =  np.random.rand(2)
distribution4 = 'VIC'
vic_inputs = [VIC_dict,vic_id_dict, fract_dict]

# Instantiate the 'LandslideProbability' component to work on this grid,
# and run it.
#Uniform
#LS_prob1 = LandslideProbability(grid,number_of_iterations=n,
#    groundwater__recharge_distribution=distribution1,
#    groundwater__recharge_min_value=Remin_value,
#    groundwater__recharge_max_value=Remax_value)
#LS_prob1.calculate_landslide_probability()
#lognormal_uniform
#LS_prob2 = LandslideProbability(grid,number_of_iterations=n,
#    groundwater__recharge_distribution=distribution2,
#    groundwater__recharge_mean=Remean,
#    groundwater__recharge_standard_deviation=Restandard_deviation)
#LS_prob2.calculate_landslide_probability()
#lognormal_spatial
#LS_prob3 = LandslideProbability(grid,number_of_iterations=n,
#    groundwater__recharge_distribution=distribution3,
#    groundwater__recharge_mean=Remean3,
#    groundwater__recharge_standard_deviation=Restandard_deviation3)
#LS_prob3.calculate_landslide_probability()
#VIC
LS_prob4 = LandslideProbability(grid,number_of_iterations=n,
    groundwater__recharge_distribution=distribution4,
    groundwater__recharge_vic_inputs=vic_inputs)
LS_prob4.calculate_landslide_probability()



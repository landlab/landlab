#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Oct 30 17:25:27 2020

@author: dylanward
"""

""" Tests based on doctests for original version """

import numpy as np
from landlab import RasterModelGrid
from landlab.components import ExponentialWeatherer, ExponentialWeathererIntegrated

mg = RasterModelGrid((5, 5))
soilz = mg.add_zeros("soil__depth", at="node")
soilrate = mg.add_ones("soil_production__rate", at="node")

expw = ExponentialWeatherer(mg)
expw.calc_soil_prod_rate()
np.allclose(mg.at_node['soil_production__rate'], 1.)
    
""" Tests for new version """

# Drop-in compatiblity
expw2 = ExponentialWeathererIntegrated(mg)
expw2.calc_soil_prod_rate()
print( mg.at_node['soil_production__rate'] )
print( np.allclose(mg.at_node['soil_production__rate'], 1.) )

expw2.run_one_step()
print('implicit soln', mg.at_node['soil_production__dt_weathered_depth'])
print( np.allclose(mg.at_node['soil_production__dt_weathered_depth'], 0.) )

# With a timestep
dt = 1000
expw2.run_one_step(dt)

print('euler soln', dt * mg.at_node['soil_production__rate'])
print('implicit soln', mg.at_node['soil_production__dt_weathered_depth'])
print( np.allclose(mg.at_node['soil_production__dt_weathered_depth'][mg.core_nodes], 6.90875478) )
print( np.allclose(mg.at_node['soil_production__dt_produced_depth'][mg.core_nodes], 6.90875478) )


# Different densities
expw3 = ExponentialWeathererIntegrated(mg,soil_production__maximum_rate=1.0, 
                                       soil_production__decay_depth=1.0,
                                       soil_production__expansion_factor=1.3)

expw3.run_one_step(dt)

print('weathered', mg.at_node['soil_production__dt_weathered_depth'])
print('produced', mg.at_node['soil_production__dt_produced_depth'])
print( np.allclose(mg.at_node['soil_production__dt_weathered_depth'][mg.core_nodes], 5.51606806) )
print( np.allclose(mg.at_node['soil_production__dt_produced_depth'][mg.core_nodes], 7.17088848) )


# Variable depths
soilz[mg.node_y > 2] =1
expw3.run_one_step(dt)

print('weathered', mg.at_node['soil_production__dt_weathered_depth'])
print('produced', mg.at_node['soil_production__dt_produced_depth'])


""" Fun with setters and getters """
print(expw2.maximum_weathering_rate, expw2._w0)
expw2.maximum_weathering_rate = 0.5


print(expw2.maximum_weathering_rate, expw2._w0)






#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Mon Jun  5 15:46:15 2017

@author: Jordan
"""

import landlab
from landlab.components import KinematicWaveRengers, DetachmentLtdErosion, SoilInfiltrationGreenAmpt
from landlab.io import read_esri_ascii
import numpy as np
from matplotlib import pyplot as plt
from landlab.plot import imshow_grid


(rmg, z) = read_esri_ascii('/Users/Jordan/Desktop/py Files/topo_as1_100cm_fr.txt', name='topographic__elevation', halo=1)
#rmg.set_watershed_boundary_condition(z, -9999.)
rmg.set_watershed_boundary_condition_outlet_id(9049, z, nodata_value = -9999.)

rmg['node']['topographic__elevation'] = z
rmg['node']['surface_water__depth'] = np.ones(rmg.number_of_nodes) * 1.e-8

kwr = KinematicWaveRengers(rmg, mannings_n = 0.02, max_courant=0.05)
dt = 1e-3

#rmg.set_watershed_boundary_condition(z)
rainrate_ms =  2.77778e-5

t = 0
maxt =60.0

stage = []
time = []

counter = 0.
samplenode = rmg.neighbors_at_node[281]
samplenode_new = samplenode[0]

while t < maxt:

    kwr.run_one_step(dt, rainfall_intensity=rainrate_ms)

    stage.append(rmg.at_node['surface_water__depth'][samplenode_new])
    time.append(t)
    
    t += dt
    dt = kwr.calc_grads_and_timesteps(False, False)


z_new = z.reshape(rmg.number_of_node_rows, rmg.number_of_node_columns)
h_new = rmg.at_node['surface_water__depth'].reshape(rmg.number_of_node_rows, rmg.number_of_node_columns)
hrow = h_new[6, :]


plt.plot(range(0, 174), hrow*1000)
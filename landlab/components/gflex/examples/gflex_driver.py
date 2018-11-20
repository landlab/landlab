# -*- coding: utf-8 -*-
"""
A driver for our version of AW's gFlex component.

Created on Fri Feb 20 11:17:52 2015

@author: danhobley
"""
from __future__ import print_function

import numpy as np
import pylab

from landlab import ModelParameterDictionary, RasterModelGrid
from landlab.components.gflex.flexure import gFlex
from landlab.plot.imshow import imshow_node_grid

inputs = ModelParameterDictionary("./AW_gflex_params.txt")
nrows = inputs.read_int("nrows")
ncols = inputs.read_int("ncols")
dx = inputs.read_float("dx")
dt = inputs.read_float("dt")
time_to_run = inputs.read_float("run_time")
init_elev = inputs.read_float("init_elev")

mg = RasterModelGrid(nrows, ncols, dx)

# create the fields in the grid
mg.add_zeros("topographic__elevation", at="node")
z = mg.zeros(at="node") + init_elev
mg["node"]["topographic__elevation"] = z + np.random.rand(len(z)) / 1000.

# make some surface load stresses in a field to test
mg.at_node["surface_load__stress"] = np.zeros(nrows * ncols, dtype=float)
square_qs = mg.at_node["surface_load__stress"].view().reshape((nrows, ncols))
square_qs[10:40, 10:40] += 1.e6

# instantiate:
gf = gFlex(mg, "./AW_gflex_params.txt")

# perform the loop:
elapsed_time = 0.  # total time in simulation
while elapsed_time < time_to_run:
    print(elapsed_time)
    if elapsed_time + dt > time_to_run:
        print("Short step!")
        dt = time_to_run - elapsed_time
    gf.flex_lithosphere()
    elapsed_time += dt

pylab.figure(1)
im = imshow_node_grid(mg, "topographic__elevation")  # display a colored image

pylab.figure(2)
im = imshow_node_grid(mg, "lithosphere_surface__elevation_increment")

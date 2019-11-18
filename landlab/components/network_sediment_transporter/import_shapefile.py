# -*- coding: utf-8 -*-
"""
Created on Wed May 22 11:48:58 2019

This code reads in data from a shapefile.

This code has been adapted from Katy's "Example import of Network Data.ipynb":
    https://github.com/landlab/tutorials/blob/barnhark/nmg_tutorials/network_model_grid/Example%20import%20of%20Network%20Data.ipynb

@author: Jon Czuba, Allison Pfeiffer, Katy Barnhart
"""
import os
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.collections import LineCollection
    
from landlab.io import read_shapefile

# %%   
DATA_DIR = './data/'


file = os.path.join(DATA_DIR, 'Test_Network.shp')
grid = read_shapefile(file)

# %%
# for plotting and testing
x_of_polylines = grid['link']['x_of_polyline']
y_of_polylines = grid['link']['y_of_polyline']

segments = []

for i in range(len(x_of_polylines)):
    x = np.array(x_of_polylines[i])
    y = np.array(y_of_polylines[i])
    segment = np.array((x, y)).T
    segments.append(segment)


from landlab.plot import graph
fig, ax = plt.subplots(figsize=(8,8), dpi=300)


graph.plot_links(grid, color='c', linestyle='solid', with_id=False,
               as_arrow=False, linewidth=1)
graph.plot_nodes(grid, color='r', with_id=False, markersize=1)

line_segments = LineCollection(segments, color='b', linewidth=0.5)
ax.add_collection(line_segments)

fig.show()

# %%
plt.figure(0)
graph.plot_graph(grid, at="node,link")

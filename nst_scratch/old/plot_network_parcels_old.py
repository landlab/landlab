# -*- coding: utf-8 -*-
"""
Created on Fri Oct 18 13:08:03 2019

@author: jczuba
"""


import matplotlib.pyplot as plt
from matplotlib.collections import LineCollection

# %% Plot a background network
# We need to set the plot limits, they will not autoscale
fig, ax = plt.subplots()
ax.set_xlim(grid.x_of_node.min(), grid.x_of_node.max())
ax.set_ylim(grid.y_of_node.min(), grid.y_of_node.max())

line_segments = LineCollection(segments, color="b", linewidth=0.5)
ax.add_collection(line_segments)

plt.show()

# %% Plot parcels at X Y location on network with color and/or size given by
# an attribute

# convert location in link to X,Y along squiggly line pathway
# -> convert_location_in_link_2_XY.py does this
# HELP HELP -- need to make this into a function or short term loop through here

# plot X,Y point on graph created above and color point according to a
# certain attribute of the parcel or the link in which the parcel resides

# HELP HELP -- We can do this below, but I can't figure out how to get it on
# the plot above. Furthermore, we can select subsets of parcels for plotting
# with different colors/sizes depending on a specific attribute(s). But first
# I need to figure out how to get it all on the same plot.

plt.plot(parcel_x, parcel_y)
plt.show()

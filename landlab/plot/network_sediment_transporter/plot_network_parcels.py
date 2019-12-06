# -*- coding: utf-8 -*-
"""
Created on Fri Oct 18 11:57:55 2019

This code plots the network and parcels and colors each
parcel according to a link or parcel attribute.

@author: jczuba
"""
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.collections import LineCollection
from landlab.plot.network_sediment_transporter import locate_parcel_xy

# this would plot lines from saved shapefile,
# then calculate distance down each grid segment,
# then place each parcel along the squiggly shapefile line correctly with
# the correct color, size, etc.
#
# have multiple options to save or return plot instance.

def plot_network_parcels(grid, parcels, parcel_time, parcel_color, parcel_size, *args, **kwargs):
   
    # This code puts the polylines into a variable segments...
    # I am not sure these next 8 lines are necessary here,
    # but I don't know how to do this differently and/or better. 
    x_of_polylines = grid['link']['x_of_polyline']
    y_of_polylines = grid['link']['y_of_polyline']

    segments = []

    for i in range(len(x_of_polylines)):
        x = np.array(x_of_polylines[i])
        y = np.array(y_of_polylines[i])
        segment = np.array((x, y)).T
        segments.append(segment)

    # Determine X Y location of each parcel
    # Initialize X Y array
    X = np.ones(len(parcels.dataset.element_id))-1
    Y = np.ones(len(parcels.dataset.element_id))-1
    
    # Locate parcel XY for each parcel at a particular time
    for parcel_number in range(len(parcels.dataset.item_id)):

        XY=locate_parcel_xy(grid, parcels, parcel_time, parcel_number)
                
        X[parcel_number]=XY[0]
        Y[parcel_number]=XY[1]
           
    # plot X,Y point on delineated network and color/size point according to a
    # certain attribute of the parcel or the link in which the parcel resides

    # we can select subsets of parcels for plotting 
    # with different colors/sizes depending on a specific attribute(s).
    
    # Plot a background network and parcels
    fig, ax = plt.subplots()
    # We need to set the plot limits, they will not autoscale
    ax.set_xlim(grid.x_of_node.min(), grid.x_of_node.max())
    ax.set_ylim(grid.y_of_node.min(), grid.y_of_node.max())
    # Plot parcels
    plt.scatter(X, Y, s=parcel_size, c=parcel_color, alpha=0.5)
    # ^HELP HELP: I can not figure out how to add a colorbar and plot the dots on top of the lines

    # Plot background network    
    line_segments = LineCollection(segments, color='b', linewidth=0.5)
    ax.add_collection(line_segments)
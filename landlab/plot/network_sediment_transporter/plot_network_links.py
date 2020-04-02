# -*- coding: utf-8 -*-
"""
Created on Fri Oct 18 11:57:55 2019

This code plots the network and colors each link according to a link attribute.

@author: jczuba
"""
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.collections import LineCollection


def plot_network_links(grid, link_attribute, *args, **kwargs):
    # Need to verify that the link_attribute is valid,
    # error messages if that's not the case

    if np.size(link_attribute) != grid.number_of_links:
        msg = "Plot Network Links: attribute must be of shape (links,)"
        raise ValueError(msg)

    # plot attribute as colored link of network

    # We need to set the plot limits, they will not autoscale
    fig, ax = plt.subplots()
    ax.set_xlim(grid.x_of_node.min(), grid.x_of_node.max())
    ax.set_ylim(grid.y_of_node.min(), grid.y_of_node.max())

    # This code puts the polylines into a variable segments...
    # I am not sure these next 8 lines are necessary here,
    # but I don't know how to do this differently and/or better.
    x_of_polylines = grid["link"]["x_of_polyline"]
    y_of_polylines = grid["link"]["y_of_polyline"]

    segments = []

    for i in range(len(x_of_polylines)):
        x = np.array(x_of_polylines[i])
        y = np.array(y_of_polylines[i])
        segment = np.array((x, y)).T
        segments.append(segment)

    line_segments = LineCollection(segments)
    line_segments.set_array(link_attribute)
    ax.add_collection(line_segments)
    fig.colorbar(line_segments)

    plt.show()

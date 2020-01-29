# -*- coding: utf-8 -*-
"""
Created on Fri Oct 18 13:08:22 2019

This code takes a link number and an attribute of the links on the grid and
(1) determines the pathway from that link to the outlet, (2) compiles the
desired attribute of the links along that pathway, and (3) creates a 2-D line
plot of the desired attributes along the pathway from the input link to the
outlet.

@author: jczuba
"""

import numpy as np
import matplotlib.pyplot as plt


def plot_pathway_values(start_link_number, link_attribute, fd, grid, *args, **kwargs):
    # Need to verify that the link_attribute is valid,
    # error messages if that's not the case

    if start_link_number < 0:
        msg = "Plot Pathway Values: link number must be >= 0"
        raise ValueError(msg)

    if start_link_number > max(grid.at_link["GridID"]):
        msg = "Plot Pathway Values: link number must be an existing link"
        raise ValueError(msg)

    # initialize variable storing links indicies along the pathway and include starting link
    link_pathway = []
    link_pathway.append(start_link_number)

    link_number = start_link_number

    downstream_link_id = fd.link_to_flow_receiving_node[
            fd.downstream_node_at_link()[link_number]
            ]

    # determine the pathway (set of links) from that link to the outlet
    while downstream_link_id != grid.BAD_INDEX:

        # determine downstream link
        downstream_link_id = fd.link_to_flow_receiving_node[
                fd.downstream_node_at_link()[link_number]
        ]

        # append downstream link id to pathway
        link_pathway.append(downstream_link_id)

        # "move" to the downstream link
        link_number = downstream_link_id

    # remove bad index value at the end of the pathway
    link_pathway.remove(-1)

    # along the pathway, index link attribute values
    # create x,y variables that create a "stairstep" appearance so each link shows that value for the entire length of its link

    x = []
    y = []
    x = np.append(x,0)
    y = np.append(y,0)

    current_distance = 0


    for i in range(len(link_pathway)):
        x = np.append(x,current_distance)
        y = np.append(y,grid.at_link[link_attribute][link_pathway[i]])
        x = np.append(x,current_distance+grid.at_link["Length"][link_pathway[i]])
        y = np.append(y,grid.at_link[link_attribute][link_pathway[i]])

        current_distance = current_distance+grid.at_link["Length"][link_pathway[i]]

    x = np.append(x,current_distance)
    y = np.append(y,0)

    # create a 2-d line plot of some specified attribute (y-axis) along that
    # pathway (distance downstream; x-axis)
    plt.figure(1)
    plt.plot(x/1000, y) #assume x is m and /1000 converts to km

    plt.xlabel("Distance downstream, km")
    plt.ylabel(link_attribute)
    plt.title(start_link_number)

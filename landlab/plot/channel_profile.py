#! /usr/bin/env python
'''
Plotting functions to extract and plot channel long profiles.
Call all three functions in sequence from the main code.

The functions will return the long profile nodes, return distances upstream of
those nodes, and plot the long profiles, respectively. The former two are, by
necessity, ragged lists of arrays - the number of channels, then the nodes in
that channel.

This module selects channels by taking the largest possible drainages crossing
the grid boundaries. You can specify how many different channels it handles
using the number_of_channels parameter in the channel_nodes function (default 
is 1). This may lead to strange outputs if the drainage structure of the output
changes mid-run (e.g., channel piracy). This may be modified in the future.

DEJH, March 2014.
Not yet tested!
'''

import numpy
try:
    import matplotlib.pyplot as plt
except ImportError:
    import warnings
    warnings.warn('matplotlib not found', ImportWarning)

def channel_nodes(grid, steepest_nodes, drainage_area, upstream_ID_order, flow_receiver, number_of_channels=1, threshold=2):
    assert threshold > 1
    boundary_nodes = grid.get_boundary_nodes()
    #top_two_pc = len(boundary_nodes)//50
    #starting_nodes = boundary_nodes[numpy.argsort(drainage_area[boundary_nodes])[-top_two_pc:]]
    starting_nodes = boundary_nodes[numpy.argsort(drainage_area[boundary_nodes])[-number_of_channels:]]
    
    profile_IDs = []
    for i in starting_nodes:
        j=i
        data_store = []
        while 1:
            data_store.append(j)
            supplying_nodes = numpy.where(flow_receiver == j)[0]
            supplying_nodes = supplying_nodes[numpy.where(supplying_nodes != i)]
            max_drainage = numpy.argmax(drainage_area[supplying_nodes])
            if drainage_area[supplying_nodes[max_drainage]] < threshold:
                break
            else:
                j = supplying_nodes[max_drainage]
        profile_IDs.append(numpy.array(data_store))
    return profile_IDs

def get_distances_upstream(grid, len_node_arrays, profile_IDs, links_to_flow_receiver):
    #defined_flow_receivers = numpy.greater_equal(links_to_flow_receiver,-1)
    #flow_link_lengths = numpy.empty(len_node_arrays, dtype=float)
    #flow_link_lengths[defined_flow_receivers] = grid.link_length[links_to_flow_receiver[defined_flow_receivers]]
    #print numpy.sum(defined_flow_receivers)
    distances_upstream = []
    for i in xrange(len(profile_IDs)):
        data_store = []
        total_distance=0.
        data_store.append(total_distance)
        for j in xrange(len(profile_IDs[i])-1):
            total_distance += grid.link_length[links_to_flow_receiver[profile_IDs[i][j+1]]]
            data_store.append(total_distance)
        distances_upstream.append(numpy.array(data_store))
    return distances_upstream

def plot_profiles(distances_upstream, profile_IDs, elevations):
    for i in xrange(len(profile_IDs)):
        the_nodes = profile_IDs[i]
        plt.plot(distances_upstream[i], elevations[the_nodes])

#! /usr/bin/env python
"""Extract and plot channel long profiles.

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
"""
# DEJH, March 2014.
from six.moves import range

import numpy
try:
    import matplotlib.pyplot as plt
except ImportError:
    import warnings
    warnings.warn('matplotlib not found', ImportWarning)


def channel_nodes(grid, steepest_nodes, drainage_area, flow_receiver, number_of_channels=1, threshold=None):
    if threshold == None:
        threshold = 2. * numpy.amin(grid.area_of_cell)
    boundary_nodes = grid.boundary_nodes
    #top_two_pc = len(boundary_nodes)//50
    #starting_nodes = boundary_nodes[numpy.argsort(drainage_area[boundary_nodes])[-top_two_pc:]]
    starting_nodes = boundary_nodes[numpy.argsort(
        drainage_area[boundary_nodes])[-number_of_channels:]]

    profile_IDs = []
    for i in starting_nodes:
        j = i
        data_store = []
        while 1:
            data_store.append(j)
            supplying_nodes = numpy.where(flow_receiver == j)[0]
            supplying_nodes = supplying_nodes[
                numpy.where(supplying_nodes != i)]
            max_drainage = numpy.argmax(drainage_area[supplying_nodes])
            if drainage_area[supplying_nodes[max_drainage]] < threshold:
                break
            else:
                j = supplying_nodes[max_drainage]
        profile_IDs.append(numpy.array(data_store))
    return profile_IDs


def get_distances_upstream(grid, len_node_arrays, profile_IDs,
                           links_to_flow_receiver):
    distances_upstream = []
    for i in range(len(profile_IDs)):
        data_store = []
        total_distance = 0.
        data_store.append(total_distance)
        for j in range(len(profile_IDs[i]) - 1):
            total_distance += grid._length_of_link_with_diagonals[
                links_to_flow_receiver[profile_IDs[i][j + 1]]]
            data_store.append(total_distance)
        distances_upstream.append(numpy.array(data_store))
    return distances_upstream


def plot_profiles(distances_upstream, profile_IDs, elevations):
    for i in range(len(profile_IDs)):
        the_nodes = profile_IDs[i]
        plt.plot(distances_upstream[i], elevations[the_nodes])


def analyze_channel_network_and_plot(grid, elevations='topographic__elevation',
                                     drainage_area='drainage_area',
                                     flow_receiver='flow__receiver_node',
                                     links_to_flow_receiver='flow__link_to_receiver_node',
                                     number_of_channels=1,
                                     starting_nodes=None,
                                     threshold=None):
    """analyze_channel_network_and_plot(grid, elevations='topographic__elevation',
                                     drainage_area='drainage_area',
                                     flow_receiver='flow__receiver_node',
                                     links_to_flow_receiver='flow__link_to_receiver_node',
                                     number_of_channels=1,
                                     starting_nodes=None,
                                     threshold=None)

    This function wraps the other three present here, and allows a single-line
    call to plot long profiles.
    As typical elsewhere, the inputs can be field names or arrays.

    Note the key new parameter starting_nodes. This (optionally) can be a
    Python list of node IDs marking the start of each profile. If it is not
    provided, the profiles with the largest terminal drainage area will be used
    instead.

    Returns a tuple, containing:
        - the list of arrays profile_IDs.
        - the list of arrays dists_upstr, the distances from the final, lowest
            node in the network.
        Both lists are number_of_channels long.
        -
    """
    internal_list = [
        0, 0, 0, 0]  # we're going to put the input arrays in here; must be a better way but this will do
    inputs = (elevations, drainage_area, flow_receiver, links_to_flow_receiver)
    for i in range(4):
        j = inputs[i]
        if type(j) == str:
            internal_list[i] = grid.at_node[j]
        else:
            assert j.size == grid.number_of_nodes, "Inputs must be field names or nnode-long numpy arrays!"
            internal_list[i] = j

    if starting_nodes == None:
        profile_IDs = channel_nodes(grid, None, internal_list[1], internal_list[
                                    2], number_of_channels, threshold)
    else:
        assert len(
            starting_nodes) == number_of_channels, "Length of starting_nodes must equal the number_of_channels!"
        if threshold == None:
            threshold = 2. * numpy.amin(grid.area_of_cell)
        profile_IDs = []
        for i in starting_nodes:
            j = i
            data_store = []
            while 1:
                data_store.append(j)
                supplying_nodes = numpy.where(flow_receiver == j)[0]
                supplying_nodes = supplying_nodes[
                    numpy.where(supplying_nodes != i)]
                max_drainage = numpy.argmax(internal_list[1][supplying_nodes])
                if internal_list[1][supplying_nodes[max_drainage]] < threshold:
                    break
                else:
                    j = supplying_nodes[max_drainage]
            profile_IDs.append(numpy.array(data_store))

    dists_upstr = get_distances_upstream(
        grid, internal_list[1].size, profile_IDs, internal_list[3])
    plot_profiles(dists_upstr, profile_IDs, internal_list[0])

    return (profile_IDs, dists_upstr)

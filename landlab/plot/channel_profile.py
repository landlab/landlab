#! /usr/bin/env python
"""Extract and plot channel long profiles.

Plotting functions to extract and plot channel long profiles.
Call all three functions in sequence from the main code.

The functions will return the long profile nodes, return distances upstream of
those nodes, and plot the long profiles, respectively. The former two are, by
necessity, ragged lists of arrays - the number of channels, then the nodes in
that channel.

You can specify how many different channels it handles
using the number_of_channels parameter in the channel_nodes function (default
is 1).

Two options exist for controlling which channels are plotted. Set 
main_channel_only = True (default) to plot only the largest drainage leading to 
model grid outlet node. main_channel_only = False will plot all channels with 
drainage below the threshold (default = 2 * grid cell area).

"""
# DEJH, March 2014.
from six.moves import range

import numpy
try:
    import matplotlib.pyplot as plt
except ImportError:
    import warnings
    warnings.warn('matplotlib not found', ImportWarning)
from landlab.plot import imshow_grid

def get_channel_segment(i, flow_receiver, drainage_area, threshold, main_channel_only):
    """ Get channel segment and return additional nodes to process. 
    
    If main_channel_only == False, the channel segment will be the largest
    drainage area segement moving upstream until drainage area is below the 
    threshold. If main_channel_only == True, the channel segment will progress
    upstream until two (or more) nodes draining into the channel have drainage
    area greater than the threshold. These nodes will be returned in 
    nodes_to_process so that they can be processed into their own channel 
    segments.
    
    Parameters
    ----------
    i: int, node idea of start of channel segment
    flow_receiver: nnode array, node id of flow recievers at node
    drainage_area: nnode array, drainage area at node
    threshold: float, threshold drainage area for the end of a channel
    main_channel_only: boolean
    
    Returns
    ----------
    channel_segment, list of nodes in the channel segment
    nodes_to_process, list of nodes to process. Empty list unless 
                main_channel_only = True
    """
    j = i
    channel_segment = []
    channel_upstream = True
    
    # add the reciever of j to the channel segment if it is not j. 
    recieving_node = flow_receiver[j]
    if recieving_node != j:
        channel_segment.append(recieving_node)
        
    while channel_upstream:
        
        # add the new node to the channel segment
        channel_segment.append(j)
        
        # get supplying nodes
        supplying_nodes = numpy.where(flow_receiver == j)[0]
        
        # remove supplying nodes that are the outlet node
        supplying_nodes = supplying_nodes[
            numpy.where(supplying_nodes != i)]
        
        # if only adding the biggest channel, continue upstream choosing the
        # largest node until no more nodes remain.      
        if main_channel_only:

            max_drainage = numpy.argmax(drainage_area[supplying_nodes])
            if drainage_area[supplying_nodes[max_drainage]] < threshold:
                nodes_to_process = []
                channel_upstream = False
            else:
                j = supplying_nodes[max_drainage]
                
        # if considering multiple channel segments, continue upstream until
        # there are two or more donors with sufficient discharge, then break, 
        # returning those nodes as starting points. 
        else:
            
            # get all upstream drainage areas
            upstream_das = drainage_area[supplying_nodes]
            
            # if no nodes upstream exceed the threshold, exit
            if numpy.sum(upstream_das > threshold) == 0:
                nodes_to_process = []
                channel_upstream = False
                
            # otherwise
            else:
                # if only one upstream node exceeds the threshold, proceed up
                # the channel. 
                if numpy.sum(upstream_das > threshold) == 1:
                    max_drainage = numpy.argmax(drainage_area[supplying_nodes])
                    j = supplying_nodes[max_drainage]
                # otherwise provide the multiple upstream nodes to be processed
                # into a new channel. 
                else:
                    nodes_to_process = supplying_nodes[upstream_das>threshold]
                    channel_upstream = False

    return (channel_segment, nodes_to_process) 

def channel_nodes(grid, starting_nodes, drainage_area, flow_receiver, number_of_channels=1, threshold=None, main_channel_only=True):
    """ 
    Create the profile_IDs data structure for channel network. 
    
    Parameters
    ----------
    grid, model grid instance. 
    starting_nodes, node ids of starting nodes
    drainage_area: nnode array, drainage area at node
    flow_receiver: nnode array, node id of flow recievers at node
    number_of_channels: int, optional, default = 1. Number of starting nodes to use
    threshold: float, optional default = 2*cell_area
        threshold drainage area for the end of a channel
    main_channel_only: boolean
    
    Returns
    ----------
    profile_IDs, channel segment datastructure. List of lists in which each 
        list is the node IDs of a channel segment. 

    """
    if threshold == None:
        threshold = 2. * numpy.amin(grid.area_of_cell)
    
    boundary_nodes = grid.boundary_nodes
    #top_two_pc = len(boundary_nodes)//50
    #starting_nodes = boundary_nodes[numpy.argsort(drainage_area[boundary_nodes])[-top_two_pc:]]
    if starting_nodes is None:
        starting_nodes = boundary_nodes[numpy.argsort(
            drainage_area[boundary_nodes])[-number_of_channels:]]

    if main_channel_only:
        profile_IDs = []
        for i in starting_nodes:
            (channel_segment, nodes_to_process) = get_channel_segment(i, flow_receiver, drainage_area, threshold, main_channel_only)
            profile_IDs.append(numpy.array(channel_segment))
    
    else:
        profile_IDs = []
        for i in starting_nodes:
            queue = [i]
            while len(queue) > 0:
                node_to_process = queue.pop(0)
                (channel_segment, nodes_to_process) = get_channel_segment(node_to_process, flow_receiver, drainage_area, threshold, main_channel_only)
                
                profile_IDs.append(numpy.array(channel_segment))
                queue.extend(nodes_to_process)

    return profile_IDs

def get_distances_upstream(grid, len_node_arrays, profile_IDs,
                           links_to_flow_receiver):
    """
    Get distances upstream for the profile_IDs datastructure.
    
    Parameters
    ----------
    grid, model grid instance. 
    len_node_arrays, int, length of node arrays
    profile_IDs: profile_IDs datastructure
    links_to_flow_receiver: nnode array, link id of the flow link to reciever
        node. 
    """
    distances_upstream = []
    end_distances = {profile_IDs[0][0]: 0}
    
    # for each profile
    for i in range(len(profile_IDs)):
        
        profile = profile_IDs[i]
        starting_node = profile[0]
        total_distance = end_distances[starting_node]

    
        data_store = []
        data_store.append(total_distance)
        
        # itterate up the profile
        for j in range(len(profile) - 1):
            total_distance += grid._length_of_link_with_diagonals[
                links_to_flow_receiver[profile[j + 1]]]
            data_store.append(total_distance)
        
        
        distances_upstream.append(numpy.array(data_store))
        end_distances[profile[-1]] = total_distance
            
    return distances_upstream


def plot_profiles(distances_upstream, profile_IDs, quanity):
    """
    Plot distance-upstream vs arbitrary quantiy
    
    Parameters
    ----------
    distances_upstream, distances upstream datastructure
    profile_IDs: profile_IDs datastructure
    quanity: nnode array, of at-node-quantity to plot against distance upstream.
    
    """
    for i in range(len(profile_IDs)):
        the_nodes = profile_IDs[i]
        plt.plot(distances_upstream[i], quanity[the_nodes])

def plot_channels_in_map_view(grid, profile_IDs, field='topographic__elevation',  **kwargs):
    """
    Plot channel locations in map view
    
    Parameters
    ----------
    grid, model grid instance. 
    field, name or nnode long array to plot with imshow_grid
    profile_IDs: profile_IDs datastructure
    **kwargs: additional parameters to pass to imshow_grid
    
    """
    imshow_grid(grid, field, **kwargs)
    for profile in profile_IDs:
        plt.plot(grid.x_of_node[profile], grid.y_of_node[profile])
        
def analyze_channel_network_and_plot(grid, elevations='topographic__elevation',
                                     drainage_area='drainage_area',
                                     flow_receiver='flow__receiver_node',
                                     links_to_flow_receiver='flow__link_to_receiver_node',
                                     number_of_channels=1,
                                     main_channel_only = True, 
                                     starting_nodes=None,
                                     threshold=None):
    """analyze_channel_network_and_plot(grid, elevations='topographic__elevation',
                                     drainage_area='drainage_area',
                                     flow_receiver='flow__receiver_node',
                                     links_to_flow_receiver='flow__link_to_receiver_node',
                                     number_of_channels=1,
                                     main_channel_only = True, 
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
    if starting_nodes is not None:
        assert len(starting_nodes) == number_of_channels, "Length of starting_nodes must equal the number_of_channels!"
    
    profile_IDs = channel_nodes(grid, starting_nodes, internal_list[1], internal_list[
                                    2], number_of_channels, threshold, main_channel_only)
    
    dists_upstr = get_distances_upstream(
        grid, internal_list[1].size, profile_IDs, internal_list[3])
    
    plot_profiles(dists_upstr, profile_IDs, internal_list[0])

    return (profile_IDs, dists_upstr)

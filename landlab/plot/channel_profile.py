#! /usr/bin/env python
"""Extract and plot channel long profiles.

Plotting functions to extract and plot channel long profiles.

Two options for use are available. For an integrated single-line call, use
analyze_channel_network_and_plot(). This function wraps the other three
functions and allows a single-line call to plot long profiles. First it uses
channel_nodes to get the nodes belonging to the channels. Then it uses
get_distances_upstream to get distances upstream. Finally it uses
plot_profiles to make a plot.


You can specify how many differents stream networks it handles using the
number_of_channels parameter in the channel_nodes function (default is 1). The
specific node ids for the beginning of the chanel can be passed with the
keyword argument starting_nodes. If it is not specified, then the stream
networks(s) will be chosed based on largest terminal drainage area.


Two options exist for controlling which channels within a given stream network
are plotted. Set main_channel_only = True (default) to plot only the largest
drainage leading that network's outlet node. main_channel_only = False will
plot all channels within each network with drainage below the threshold value
(default = 2 * grid cell area).

The functions will return the profile datastructure profile_structure and the
distance upstream datastructure distances_upstream. rofile structure is a list
of length number_of_channels. Each element of profile_structure is itself a
list of length number of stream  segments that drain to each of the starting
nodes. Each stream segment list contains the node ids of a stream segment from
downstream to upstream. distances_upstream provides the equivalent structure
but provides the distance upstream rather than the node id.

Examples
---------

    Start by importing necessary modules

    >>> import numpy as np
    >>> np.random.seed(42)
    >>> from landlab import RasterModelGrid
    >>> from landlab.components import FlowAccumulator, FastscapeEroder
    >>> from landlab.plot import analyze_channel_network_and_plot

    Construct a grid and evolve some topography

    >>> mg = RasterModelGrid(40, 60)
    >>> z = mg.add_zeros('topographic__elevation', at='node')
    >>> z += 200 + mg.x_of_node + mg.y_of_node
    >>> mg.set_closed_boundaries_at_grid_edges(bottom_is_closed=True, left_is_closed=True, right_is_closed=True, top_is_closed=True)
    >>> mg.set_watershed_boundary_condition_outlet_id(0, z, -9999)
    >>> fa = FlowAccumulator(mg, flow_director='D8')
    >>> sp = FastscapeEroder(mg, K_sp=.0001, m_sp=.5, n_sp=1)

    Now run a simple landscape evolution model to develop a channel network.

    >>> dt = 100
    >>> for i in range(200):
    ...     fa.run_one_step()
    ...     sp.run_one_step(dt=dt)
    ...     mg.at_node['topographic__elevation'][0] -= 0.001

    Now call analyze_channel_network_and_plot in order to plot the network. Note
    in general, we'd leave create_plot in its default value of create_plot=True,
    but we can't plot in the docstring

    >>> profile_structure, distances_upstream = analyze_channel_network_and_plot(mg,
    ...                                                                          create_plot=False)

    This creates the profile structure and distances upstream structure. Since we
    used the default values this will make only one profile, the biggest channel
    in the biggest stream network in the catchemnt.

    profile_structure will be a length 1 array and will contain one array that
    is the length in number of nodes of the single longest channel

    >>> len(profile_structure) == 1
    True

    >>> len(profile_structure[0][0])
    58

    We can change the default values to get additional channels or to plot all
    of the channels in a network.

    For the next example, lets use a hexagonal grid.

    >>> from landlab import HexModelGrid
    >>> from landlab.components import DepressionFinderAndRouter, LinearDiffuser
    >>> mg = HexModelGrid(40, 20)
    >>> z = mg.add_zeros('topographic__elevation', at='node')
    >>> z += 200 + mg.x_of_node + mg.y_of_node + np.random.randn(mg.size('node'))
    >>> fa = FlowAccumulator(mg, depression_finder=DepressionFinderAndRouter)
    >>> sp = FastscapeEroder(mg, K_sp=.0001, m_sp=.5, n_sp=1)
    >>> ld = LinearDiffuser(mg, linear_diffusivity=0.0001)

    Now run a simple landscape evolution model to develop a channel network.

    >>> dt = 100
    >>> for i in range(200):
    ...     fa.run_one_step()
    ...     flooded = np.where(fa.depression_finder.flood_status==3)[0]
    ...     sp.run_one_step(dt=dt,  flooded_nodes=flooded)
    ...     ld.run_one_step(dt=dt)
    ...     mg.at_node['topographic__elevation'][0] -= 0.001

    This time we will use some non-default values. Providing the threshold in
    units of area indicates where the channel network will end.
    main_channel_only = False indicates that all parts of the stream network
    >>> profile_structure, distances_upstream = analyze_channel_network_and_plot(mg,
    ...                                                                          threshold = 100,
    ...                                                                          main_channel_only=False,
    ...                                                                          number_of_channels=3)

    This will create four channel networks. The datastructures profile_structure
    and distances_upstream will both be length 3

    >>> len(profile_structure) == len(distances_upstream) == 3
    True

    We can also specify exactly which node is the outlet for the channel.

    >>> profile_structure, distances_upstream = analyze_channel_network_and_plot(mg,
    ...                                                                          threshold = 100,
    ...                                                                          starting_nodes = [0],
    ...                                                                          number_of_channels=1)

    It is important that the length of starting nodes is the same as the value
    of number_of_channels. If this is not the case, then an error will occur.


"""


from six.moves import range

import numpy
try:
    import matplotlib.pyplot as plt
except ImportError:
    import warnings
    warnings.warn('matplotlib not found', ImportWarning)
from landlab.plot import imshow_grid
from landlab.utils.decorators import use_field_name_or_array
from landlab import RasterModelGrid

@use_field_name_or_array('node')
def _return_surface(grid, surface):
    """
    Private function to return grid fields.

    This function exists to take advantange of the 'use_field_name_or_array
    decorator which permits providing the surface as a field name or array.
    """
    return surface

def _get_channel_segment(i, flow_receiver, drainage_area, threshold, main_channel_only):
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
    i : int, required
        Node id of start of channel segment.
    flow_receiver : nnode array, required
        Node id  array of flow recievers at node.
    drainage_area : nnode array, required
        Node id array of drainage area at node
    threshold : float, required
        Threshold drainage area for determining the end of a channel.
    main_channel_only : boolean, required
        Flag that indicates if only the main channel is being calculated or if
        all channel segments with drainage area below threshold are being
        identified.

    Returns
    ----------
    channel_segment : list
        Node IDs of the nodes in the current channel segment.
    nodes_to_process, list
        List of nodes to add to the processing queue. These nodes are those
        that drain to the upper end of this channel segment. If
        main_channel_only = False this will be an empty list.
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
    grid : Landlab Model Grid instance, required
    starting_nodes : length number_of_channels itterable, optional
        Length number_of_channels itterable containing the node IDs of nodes
        to start the channel profiles from. If not provided, the default is the
        number_of_channels node IDs on the model grid boundary with the largest
        terminal drainage area
    drainage_area : nnode array, required
        Node id array of drainage area at node
    flow_receiver : nnode array, required
        Node id  array of flow recievers at node.
    number_of_channels : int, optional
        Total number of channels to plot. Default value is 1. If value is
        greater than 1 and starting_nodes is not specified, then the
        number_of_channels largest channels based on terminal drainage area
        will be used.
    threshold : float, optional
        Value to use for the minimum drainage area associated with a plotted
        channel segment. Default values is 2.0 x minimum grid cell area.
    main_channel_only : Boolean, optional
        Flag to determine if only the main channel should be plotted, or if all
        stream segments with drainage area less than threshold should be
        plotted. Default value is True.

    Returns
    ----------
    profile_structure, the channel segment datastructure.
        profile structure is a list of length number_of_channels. Each element
        of profile_structure is itself a list of length number of stream
        segments that drain to each of the starting nodes. Each stream segment
        list contains the node ids of a stream segment from downstream to
        upstream.

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
        profile_structure = []
        channel_network = []
        for i in starting_nodes:
            (channel_segment, nodes_to_process) = _get_channel_segment(i, flow_receiver, drainage_area, threshold, main_channel_only)
            channel_network.append(numpy.array(channel_segment))
            profile_structure.append(channel_network)

    else:
        profile_structure = []
        for i in starting_nodes:
            channel_network = []
            queue = [i]
            while len(queue) > 0:
                node_to_process = queue.pop(0)
                (channel_segment, nodes_to_process) = _get_channel_segment(node_to_process, flow_receiver, drainage_area, threshold, main_channel_only)
                channel_network.append(numpy.array(channel_segment))
                queue.extend(nodes_to_process)
            profile_structure.append(channel_network)

    return profile_structure

def get_distances_upstream(grid,
                           profile_structure,
                           links_to_flow_receiver):
    """
    Get distances upstream for the profile_IDs datastructure.

    Parameters
    ----------
    grid :  model grid instance, required
    length of node arrays
    profile_structure: profile_structure datastructure
    links_to_flow_receiver: nnode array, link id of the flow link to reciever
        node.

    Returns
    ----------
    distances_upstream, datastruture that mirrors profile IDs but provides the
        distance upstream.
    """
    end_distances = {}

    # set the starting values for the beginnings of each netwrok.

    for network in profile_structure:
        starting_node = network[0][0]
        end_distances[starting_node] = 0

    distances_upstream = []
    # for each network
    for network in profile_structure:

        network_values = []
        # for each segment in the network.
        for segment in network:
            starting_node = segment[0]

            total_distance = end_distances[starting_node]

            profile_values = []
            profile_values.append(total_distance)

            # itterate up the profile
            for j in range(len(segment) - 1):
                if isinstance(grid, RasterModelGrid):
                    total_distance += grid.length_of_d8[
                        links_to_flow_receiver[segment[j + 1]]]
                    profile_values.append(total_distance)
                else:
                    total_distance += grid.length_of_link[
                        links_to_flow_receiver[segment[j + 1]]]
                    profile_values.append(total_distance)
            network_values.append(numpy.array(profile_values))
            end_distances[segment[-1]] = total_distance

        distances_upstream.append(network_values)
    return distances_upstream


def plot_profiles(distances_upstream, profile_structure, quantity):
    """
    Plot distance-upstream vs arbitrary quantity, default when calling through
    analyze_channel_network_and_plot is topographic_elevation.

    Parameters
    ----------
    distances_upstream : datastructure, required
        The distances upstream datastructure created by get_distances_upstream
    profile_structure : datastructure, required
        profile_structure datastructure created by channel_nodes
    quantity : nnode array, required
        Array of  the at-node-quantity to plot against distance upstream.
    """
    # for each stream network
    for i in range(len(profile_structure)):
        network_nodes = profile_structure[i]
        network_distance = distances_upstream[i]

        # for each stream segment in the network
        for j in range(len(network_nodes)):

            # identify the nodes and distances upstream for this channel segment
            the_nodes = network_nodes[j]
            the_distances = network_distance[j]
            plt.plot(the_distances, quantity[the_nodes])


def plot_channels_in_map_view(grid, profile_structure, field='topographic__elevation',  **kwargs):
    """
    Plot channel locations in map view on a frame.

    Parameters
    ----------
    grid, model grid instance.
    field, name or nnode long array to plot with imshow_grid
    profile_IDs: profile_IDs datastructure
    **kwargs: additional parameters to pass to imshow_grid
    """
    # make imshow_grid background
    imshow_grid(grid, field, **kwargs)

    # for each stream network
    for i in range(len(profile_structure)):
        network_nodes = profile_structure[i]

        # for each stream segment in the network
        for j in range(len(network_nodes)):

            # identify the nodes and distances upstream for this channel segment
            the_nodes = network_nodes[j]
            plt.plot(grid.x_of_node[the_nodes], grid.y_of_node[the_nodes])


def analyze_channel_network_and_plot(grid,
                                     field='topographic__elevation',
                                     drainage_area='drainage_area',
                                     flow_receiver='flow__receiver_node',
                                     links_to_flow_receiver='flow__link_to_receiver_node',
                                     number_of_channels=1,
                                     main_channel_only = True,
                                     starting_nodes=None,
                                     threshold=None,
                                     create_plot=True):
    """
    Main function to analyse the channel network and make an distance upstream
    vs the quantity stored at the model grid field give by the keyword argument
    `field`.

    This function wraps the other three present here, and allows a single-line
    call to plot long profiles. First it uses channel_nodes to get the nodes
    belonging to the channels. Then it uses get_distances_upstream to get
    distances upstream. Finally it uses plot_profiles to make a plot.

    Parameters
    ----------
    grid : Landlab Model Grid instance, required
    field : string or length nnode array, optional
        Field name or array of the quantity to plot against distance upstream.
        Default value is 'topographic__elevation'.
    drainage_area : string or length nnode array, optional
        Field name or array of the drainage area of the model grid. This field
        is used to identify the boundary nodes with the largest terminal
        drainage area and to identify the end of channels using the threshold
        parameter. Default value is 'drainage_area' which will be created by
        Landlab's FlowAccumulator or FlowRouter.
    flow_receiver : string or length nnode array, optional
        Field name or array of the flow_links to reciever node of the model
        grid. Default value is 'flow__receiver_node' which will be created by
        Landlab's FlowAccumulator or FlowRouter.
    links_to_flow_receiver : string or length nnode array, optional
        Field name or array of the flow_links to reciever node of the model
        grid. Default value is 'flow__link_to_receiver_node' which will be
        created by Landlab's FlowAccumulator or FlowRouter.
    number_of_channels : int, optional
        Total number of channels to plot. Default value is 1. If value is
        greater than 1 and starting_nodes is not specified, then the
        number_of_channels largest channels based on terminal drainage area
        will be used.
    main_channel_only : Boolean, optional
        Flag to determine if only the main channel should be plotted, or if all
        stream segments with drainage area less than threshold should be
        plotted. Default value is True.
    starting_nodes : length number_of_channels itterable, optional
        Length number_of_channels itterable containing the node IDs of nodes
        to start the channel profiles from. If not provided, the default is the
        number_of_channels node IDs on the model grid boundary with the largest
        terminal drainage area
    threshold : float, optional
        Value to use for the minimum drainage area associated with a plotted
        channel segment. Default values is 2.0 x minimum grid cell area.
    create_plot : boolean, optional
        Flag to indicate if a distance-upstream vs plotted quantity plot should
        be created. Default is True.

    Returns
    ----------
    tuple, containing:
        profile_structure, the channel segment datastructure.
            profile structure is a list of length number_of_channels. Each
            element of profile_structure is itself a list of length number of
            stream segments that drain to each of the starting nodes. Each
            stream segment list contains the node ids of a stream segment from
            downstream to upstream.
        distances_upstream, the channel segment datastructure.
            A datastructure that parallels profile_structure but holds
            distances upstream instead of node IDs.
        Both lists are number_of_channels long.
    """
    # get values for the four at-node arrays needed
    field_values = _return_surface(grid, field)
    drainage_area_values = _return_surface(grid, drainage_area)
    flow_receiver_values = _return_surface(grid, flow_receiver)
    links_to_flow_receiver_values = _return_surface(grid, links_to_flow_receiver)

    # verify that the number of starting nodes is the specified number of channels
    if starting_nodes is not None:
        assert len(starting_nodes) == number_of_channels, "Length of starting_nodes must equal the number_of_channels!"

    # get the profile IDs datastructure
    profile_structure = channel_nodes(grid,
                                      starting_nodes,
                                      drainage_area_values,
                                      flow_receiver_values,
                                      number_of_channels,
                                      threshold,
                                      main_channel_only)

    # use the profile IDs to get the distances upstream
    distances_upstream = get_distances_upstream(grid,
                                                profile_structure,
                                                links_to_flow_receiver_values)

    # if requested, create plot
    if create_plot:
        plot_profiles(distances_upstream, profile_structure, field_values)

    # return the profile IDS and the disances upstream.
    return (profile_structure, distances_upstream)

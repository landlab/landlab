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

Examples
---------
>>> import numpy as np
>>> from landlab import RasterModelGrid
>>> from landlab.components.flow_routing import FlowRouter
>>> from landlab.components import DepressionFinderAndRouter
>>> from landlab.components import Space
>>> from landlab.components import FastscapeEroder
>>> np.random.seed(seed = 5000)

Define grid and initial topography:

*  5x5 grid with baselevel in the lower left corner
*  All other boundary nodes closed
*  Initial topography is plane tilted up to the upper right with
   noise

>>> mg = RasterModelGrid((5, 5), spacing=10.0)
>>> _ = mg.add_zeros('topographic__elevation', at='node')
>>> mg.at_node['topographic__elevation'] += (mg.node_y / 10. +
...     mg.node_x / 10. + np.random.rand(len(mg.node_y)) / 10.)
>>> mg.set_closed_boundaries_at_grid_edges(bottom_is_closed=True,
...                                        left_is_closed=True,
...                                        right_is_closed=True,
...                                        top_is_closed=True)
>>> mg.set_watershed_boundary_condition_outlet_id(
...     0, mg.at_node['topographic__elevation'], -9999.)
>>> fsc_dt = 100.
>>> space_dt = 100.

Instantiate Fastscape eroder, flow router, and depression finder

>>> fsc = FastscapeEroder(mg, K_sp=.001, m_sp=.5, n_sp=1)
>>> fr = FlowRouter(mg)
>>> df = DepressionFinderAndRouter(mg)

Burn in an initial drainage network using the Fastscape eroder:

>>> for x in range(100):
...     fr.run_one_step()
...     df.map_depressions()
...     flooded = np.where(df.flood_status == 3)[0]
...     fsc.run_one_step(dt=fsc_dt, flooded_nodes=flooded)
...     mg.at_node['topographic__elevation'][0] -= 0.001 # Uplift

analyze_channel_network_and_plot(grid, elevations='topographic__elevation',
                                    drainage_area='drainage_area',
                                    flow_receiver='flow__receiver_node',
                                    links_to_flow_receiver='flow__link_to_receiver_node',
                                    number_of_channels=1,
                                    main_channel_only = True,
                                    starting_nodes=None,
                                    threshold=None,
                                    create_plot=True)
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
            (channel_segment, nodes_to_process) = _get_channel_segment(i, flow_receiver, drainage_area, threshold, main_channel_only)
            profile_IDs.append(numpy.array(channel_segment))

    else:
        profile_IDs = []
        for i in starting_nodes:
            queue = [i]
            while len(queue) > 0:
                node_to_process = queue.pop(0)
                (channel_segment, nodes_to_process) = _get_channel_segment(node_to_process, flow_receiver, drainage_area, threshold, main_channel_only)

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

    Returns
    ----------
    distances_upstream, datastruture that mirrors profile IDs but provides the
        distance upstream.
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
    Plot channel locations in map view on a frame.

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
    grid : Raster Model Grid, required
    field : string or length nnode array, optional
        Field name or array of the quantity to plot against distance upstream.
        Default value is 'topographic__elevation'.
    drainage_area : string or length nnode array, optional
        Field name or array of the drainage area of the model grid.
        Default value is 'drainage_area' which will be created by Landlab's
        FlowAccumulator or FlowRouter.
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
        number_of_channels largest channels based on drainage area will be used.
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
        - profile_IDs datastructure
        - the list of arrays dists_upstr, the distances from the final, lowest
            node in the network.
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
    profile_IDs = channel_nodes(grid,
                                starting_nodes,
                                drainage_area_values,
                                flow_receiver_values,
                                number_of_channels,
                                threshold,
                                main_channel_only)

    # use the profile IDs to get the distances upstream
    dists_upstr = get_distances_upstream(grid,
                                         drainage_area_values.size,
                                         profile_IDs,
                                         links_to_flow_receiver_values)

    # if requested, create plot
    if create_plot:
        plot_profiles(dists_upstr, profile_IDs, field_values)

    # return the profile IDS and the disances upstream.
    return (profile_IDs, dists_upstr)

# coding: utf8
#! /usr/env/python
"""
"""

from six.moves import range
import numpy as np

from landlab import RasterModelGrid
from landlab.components.profiler.base_profiler import _Profiler


class ChannelProfiler(_Profiler):
    """Profile channels of a drainage network.

    It is expected that the following at-node grid fields will be present.

        'drainage_area'
        'flow__receiver_node'
        'flow__link_to_receiver_node'

    Extract and plot channel long profiles.

    Two options for use are available. For an integrated single-line call, use
    analyze_channel_network_and_plot(). This function wraps the other three
    functions and allows a single-line call to plot long profiles. First it uses
    channel_nodes to get the nodes belonging to the channels. Then it uses
    get_distances_upstream to get distances upstream. Finally it uses
    plot_profiles to make a plot.


    You can specify how many differents stream networks it handles using the
    number_of_watersheds parameter in the channel_nodes function (default is 1). The
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
    of length number_of_watersheds. Each element of profile_structure is itself a
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
    >>> from landlab.components import (FlowAccumulator,
    ...                                 FastscapeEroder,
    ...                                 ChannelProfiler)

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

    >>> profiler = ChannelProfiler(mg)
    >>> profiler.run_one_step()

    This creates the profile structure and distances upstream structure. Since we
    used the default values this will make only one profile, the biggest channel
    in the biggest stream network in the catchemnt.

    profile_structure will be a length 1 array and will contain one array that
    is the length in number of nodes of the single longest channel

    >>> len(profiler.profile_structure) == 1
    True

    >>> len(profiler.profile_structure[0][0])
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

    >>> profiler = ChannelProfiler(mg,
    ...                            threshold = 100,
    ...                            main_channel_only=False,
    ...                            number_of_watersheds=3)

    This will create four channel networks. The datastructures profile_structure
    and distances_upstream will both be length 3

    >>> len(profiler.profile_structure) == len(profiler.distances_upstream) == 3
    True

    We can also specify exactly which node is the outlet for the channel.

    >>> profiler = ChannelProfiler(mg,
    ...                            threshold = 100,
    ...                            starting_nodes = [0],
    ...                            number_of_watersheds=1)

    It is important that the length of starting nodes is the same as the value
    of number_of_watersheds. If this is not the case, then an error will occur.

    """
    def __init__(self, grid,
                       number_of_watersheds=1,
                       main_channel_only = True,
                       starting_nodes=None,
                       threshold=None):
        """Parameters
        ----------
        grid : Landlab Model Grid instance, required
        number_of_watersheds : int, optional
            Total number of watersheds to plot. Default value is 1. If value is
            greater than 1 and starting_nodes is not specified, then the
            number_of_watersheds largest watersheds based on the drainage area
            at the model grid boundary.
        main_channel_only : Boolean, optional
            Flag to determine if only the main channel should be plotted, or if all
            stream segments with drainage area less than threshold should be
            plotted. Default value is True.
        starting_nodes : length number_of_watersheds itterable, optional
            Length number_of_watersheds itterable containing the node IDs of nodes
            to start the channel profiles from. If not provided, the default is the
            number_of_watersheds node IDs on the model grid boundary with the largest
            terminal drainage area
        threshold : float, optional
            Value to use for the minimum drainage area associated with a plotted
            channel segment. Default values is 2.0 x minimum grid cell area.
        """
        super(ChannelProfiler, self).__init__(grid)

        self._distances_upstream = []
        self._main_channel_only = main_channel_only

        if threshold is None:
            threshold = 2. * np.amin(grid.area_of_cell)
        self.threshold = threshold

        # verify that the number of starting nodes is the specified number of channels
        if starting_nodes is not None:
            if len(starting_nodes) is not number_of_watersheds:
                msg = "Length of starting_nodes must equal the number_of_watersheds!"
                raise ValueError(msg)
        else:
            starting_nodes = grid.boundary_nodes[np.argsort(
                self._drainage_area[grid.boundary_nodes])[-number_of_watersheds:]]

        self._starting_nodes = starting_nodes

    @property
    def profile_structure():
        """
        profile_structure, the channel segment datastructure.
                profile structure is a list of length number_of_watersheds. Each
                element of profile_structure is itself a list of length number of
                stream segments that drain to each of the starting nodes. Each
                stream segment list contains the node ids of a stream segment from
                downstream to upstream.
        """
        return self._profile_structure
    @property
    def distances_upstream():
        """
        distances_upstream, the channel segment datastructure.
                A datastructure that parallels profile_structure but holds
                distances upstream instead of node IDs.

                Both lists are number_of_watersheds long.
        """
        return self._distances_upstream

    def run_one_step(self):
        """Calculate the channel profile datastructure and distances upstream.


        """
        # calculate the profile IDs datastructure
        self._create_profile_structure()

        # calculate the distance along profile datastructure
        self._calculate_distances_upstream()

    def _get_channel_segment(self, i):
        """Get channel segment and return additional nodes to process.

        Parameters
        ----------
        i : int, required
            Node id of start of channel segment.

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
        recieving_node = self._flow_receiver[j]
        if recieving_node != j:
            channel_segment.append(recieving_node)

        while channel_upstream:

            # add the new node to the channel segment
            channel_segment.append(j)

            # get supplying nodes
            supplying_nodes = np.where(self._flow_receiver == j)[0]

            # remove supplying nodes that are the outlet node
            supplying_nodes = supplying_nodes[
                np.where(supplying_nodes != i)]

            # if only adding the biggest channel, continue upstream choosing the
            # largest node until no more nodes remain.
            if self._main_channel_only:
                max_drainage = np.argmax(self._drainage_area[supplying_nodes])
                if self._drainage_area[supplying_nodes[max_drainage]] < self.threshold:
                    nodes_to_process = []
                    channel_upstream = False
                else:
                    j = supplying_nodes[max_drainage]

            # if considering multiple channel segments, continue upstream until
            # there are two or more donors with sufficient discharge, then break,
            # returning those nodes as starting points.
            else:

                # get all upstream drainage areas
                upstream_das = self._drainage_area[supplying_nodes]

                # if no nodes upstream exceed the threshold, exit
                if np.sum(upstream_das > self.threshold) == 0:
                    nodes_to_process = []
                    channel_upstream = False

                # otherwise
                else:
                    # if only one upstream node exceeds the threshold, proceed up
                    # the channel.
                    if np.sum(upstream_das > self.threshold) == 1:
                        max_drainage = np.argmax(self._drainage_area[supplying_nodes])
                        j = supplying_nodes[max_drainage]
                    # otherwise provide the multiple upstream nodes to be processed
                    # into a new channel.
                    else:
                        nodes_to_process = supplying_nodes[upstream_das>self.threshold]
                        channel_upstream = False

        return (channel_segment, nodes_to_process)

    def _create_profile_structure(self):
        """
        Create the profile_IDs data structure for channel network.

        The bound attribute self._profile structure is the channel segment
        datastructure. profile structure is a list of length
        number_of_watersheds. Each element of profile_structure is itself a list
        of length number of stream segments that drain to each of the starting
        nodes. Each stream segment list contains the node ids of a stream
        segment from downstream to upstream.
        """
        self._profile_structure = []

        if self._main_channel_only:
            channel_network = []
            for i in self._starting_nodes:
                (channel_segment, nodes_to_process) = self._get_channel_segment(i)
                channel_network.append(np.array(channel_segment))
                self._profile_structure.append(channel_network)

        else:
            for i in self._starting_nodes:
                channel_network = []
                queue = [i]
                while len(queue) > 0:
                    node_to_process = queue.pop(0)
                    (channel_segment, nodes_to_process) = self._get_channel_segment(node_to_process)
                    channel_network.append(np.array(channel_segment))
                    queue.extend(nodes_to_process)
                self._profile_structure.append(channel_network)

    def _calculate_distances_upstream(self):
        """
        Get distances upstream for the profile_IDs datastructure.
        """
        end_distances = {}

        # set the starting values for the beginnings of each netwrok.
        for network in self._profile_structure:
            starting_node = network[0][0]
            end_distances[starting_node] = 0

        # for each network
        for network in self._profile_structure:

            network_values = []
            # for each segment in the network.
            for segment in network:
                starting_node = segment[0]

                total_distance = end_distances[starting_node]

                profile_values = []
                profile_values.append(total_distance)

                # itterate up the profile
                for j in range(len(segment) - 1):
                    if isinstance(self._grid, RasterModelGrid):
                        total_distance += self._grid.length_of_d8[
                            self._link_to_flow_receiver[segment[j + 1]]]
                        profile_values.append(total_distance)
                    else:
                        total_distance += self._grid.length_of_link[
                            self._link_to_flow_receiver[segment[j + 1]]]
                        profile_values.append(total_distance)
                network_values.append(np.array(profile_values))
                end_distances[segment[-1]] = total_distance
            self._distances_upstream.append(network_values)

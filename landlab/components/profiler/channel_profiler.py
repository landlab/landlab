# coding: utf8
# ! /usr/env/python
"""channel_profiler.py component to create channel profiles."""

import numpy as np
from six.moves import range

from landlab import RasterModelGrid
from landlab.components.profiler.base_profiler import _BaseProfiler


class ChannelProfiler(_BaseProfiler):
    """Extract and plot the long profiles within a drainage network.

    The ChannelProfiler can work on grids that have used route-to-one or
    route-to-multiple flow directing. This component expects that all of the
    following at-node grid fields will be present:

        'drainage_area'
        'flow__receiver_node'
        'flow__link_to_receiver_node'

    These are typically created by using the FlowAccumulator component.

    Four parameters control the behavior of this component. By default, the
    component will extract the single largest channel in the drainage network
    up to a drainage area below a threshold of two times the cell area.

    Providing a value for the parameter ``threshold`` changes this value.
    Setting ``main_channel_only`` to ``False`` will plot all of the channels in
    each watershed that have drainage area below the threshold. The number of
    watersheds can be changed from one based on the parameter
    ``number_of_watersheds``. By default the ChannelProfiler will identify the
    largest watersheds and use them. If instead you want to extract the channel
    network draining to specific nodes, you can provide values for the
    parameter ``starting_nodes``. The length of ``starting_nodes`` must be
    consistent with the value of ``number_of_watersheds``.  The node IDs of the
    ``starting_nodes`` do not need to be on the domain boundaries or the
    outlets of watersheds. They do need to have drainage area above the
    indicated threshold.

    The ``run_one_step`` method extracts the channel network and stores it in
    the ``profile_structure`` and ``distances_upstream`` bound properties. To
    make a plot of distance-upstream vs an at-node value use the method
    ``plot_profiles``. To make a plot of the location of profiles in map view
    use ``plot_profiles_in_map_view``.

    Examples
    ---------

    Start by importing necessary modules

    >>> import numpy as np
    >>> np.random.seed(42)
    >>> from landlab import RasterModelGrid
    >>> from landlab.components import (FlowAccumulator,
    ...                                 FastscapeEroder,
    ...                                 ChannelProfiler)

    Construct a grid and set up components necessary to evolve some topography

    >>> mg = RasterModelGrid(40, 60)
    >>> z = mg.add_zeros('topographic__elevation', at='node')
    >>> z += 200 + mg.x_of_node + mg.y_of_node
    >>> mg.set_closed_boundaries_at_grid_edges(bottom_is_closed=True,
    ...                                        left_is_closed=True,
    ...                                        right_is_closed=True,
    ...                                        top_is_closed=True)
    >>> mg.set_watershed_boundary_condition_outlet_id(0, z, -9999)
    >>> fa = FlowAccumulator(mg, flow_director='D8')
    >>> sp = FastscapeEroder(mg, K_sp=.0001, m_sp=.5, n_sp=1)

    Now run our simple landscape evolution model to develop a channel network.

    >>> dt = 100
    >>> for i in range(200):
    ...     fa.run_one_step()
    ...     sp.run_one_step(dt=dt)
    ...     mg.at_node['topographic__elevation'][0] -= 0.001

    Next we construct the ChannelProfiler component and run ``run_one_step`` to
    extract the channel network.

    >>> profiler = ChannelProfiler(mg)
    >>> profiler.run_one_step()

    This creates the profile structure and distances upstream structure. Since
    we used the default values this will make only one profile, the biggest
    channel in the biggest stream network in the catchemnt.

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
    >>> from landlab.components import (
    ...     DepressionFinderAndRouter,
    ...     LinearDiffuser)
    >>> mg = HexModelGrid(40, 20)
    >>> z = mg.add_zeros('topographic__elevation', at='node')
    >>> z += (200
    ...       + mg.x_of_node
    ...       + mg.y_of_node
    ...       + np.random.randn(mg.size('node')))
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
    >>> profiler.run_one_step()

    This will create three channel networks. The attributes profile_structure
    and distances_upstream will both be length 3

    >>> len(profiler.profile_structure) == len(profiler.distances_upstream) == 3
    True

    We can also specify exactly which node is the outlet for the channel.

    >>> profiler = ChannelProfiler(mg,
    ...                            threshold = 300,
    ...                            starting_nodes = [0],
    ...                            number_of_watersheds=1)
    >>> profiler.run_one_step()
    >>> np.testing.assert_array_equal(
    ...     profiler.profile_structure[0][0],
    ...     np.array([ 0, 21, 22, 23, 24], dtype=np.int64))

    It is important that the length of starting nodes is the same as the value
    of number_of_watersheds. If this is not the case, then an error will occur.

    """

    _name = "ChannelProfiler"

    _input_var_names = (
        "topographic__elevation",
        "drainage_area",
        "flow__receiver_node",
        "flow__link_to_receiver_node",
    )

    _output_var_names = ()

    _var_units = {
        "topographic__elevation": "m",
        "flow__receiver_node": "-",
        "drainage_area": "m**2",
        "flow__link_to_receiver_node": "-",
    }

    _var_mapping = {
        "topographic__elevation": "node",
        "flow__receiver_node": "node",
        "drainage_area": "node",
        "flow__link_to_receiver_node": "node",
    }
    _var_doc = {
        "topographic__elevation": "Land surface topographic elevation",
        "flow__receiver_node": "Node array of receivers (node that receives flow from current node)",
        "drainage_area": "Upstream accumulated surface area contributing to the node's discharge",
        "flow__link_to_receiver_node": "Node array containing ID of link that leads from each node to its receiver, or BAD_INDEX_VALUE if no link",
    }

    def __init__(
        self,
        grid,
        stopping_field="drainage_area",
        number_of_watersheds=1,
        outlet_threshold=1e6,
        main_channel_only=True,
        starting_nodes=None,
        threshold=None,
    ):
        """Parameters
        ----------
        grid : Landlab Model Grid instance, required
        stopping_field : field name as string
            Field name to indicate a basis for stopping the channel profile.
            Default is "drainage_area".
        number_of_watersheds : int, optional
            Total number of watersheds to plot. Default value is 1. If value is
            greater than 1 and starting_nodes is not specified, then the
            number_of_watersheds largest watersheds based on the drainage area
            at the model grid boundary. If given as None, then all grid cells
            on the domain boundary with a stopping field (typically drainage
            area) greater than the outlet_threshold in area are used.
        outlet_threshold : float, optional
            Threshold for defining an outlet. Default is 1e6.
        main_channel_only : Boolean, optional
            Flag to determine if only the main channel should be plotted, or if
            all stream segments with drainage area less than threshold should
            be plotted. Default value is True.
        starting_nodes : length number_of_watersheds itterable, optional
            Length number_of_watersheds itterable containing the node IDs of
            nodes to start the channel profiles from. If not provided, the
            default is the number_of_watersheds node IDs on the model grid
            boundary with the largest terminal drainage area.
        threshold : float, optional
            Value to use for the minimum drainage area associated with a
            plotted channel segment. Default values is 2.0 x minimum grid cell
            area.
        """
        super(ChannelProfiler, self).__init__(grid, stopping_field)

        self._grid = grid

        if stopping_field in grid.at_node:
            self._stopping_field = grid.at_node[stopping_field]
        else:
            msg = (
                "a field to stop based on is a required field to run a ChannelProfiler."
            )
            raise ValueError(msg)

        if "flow__receiver_node" in grid.at_node:
            self._flow_receiver = grid.at_node["flow__receiver_node"]
        else:
            msg = "flow__receiver_node is a required field to run a ChannelProfiler."
            raise ValueError(msg)

        if "flow__link_to_receiver_node" in grid.at_node:
            self._link_to_flow_receiver = grid.at_node["flow__link_to_receiver_node"]
        else:
            msg = "flow__link_to_receiver_node is a required field to run a ChannelProfiler."
            raise ValueError(msg)

        self._main_channel_only = main_channel_only

        if threshold is None:
            threshold = 2.0 * np.amin(grid.area_of_cell)
        self.threshold = threshold

        # verify that the number of starting nodes is the specified number of channels
        if starting_nodes is not None:
            if (number_of_watersheds is not None) and (
                len(starting_nodes) is not number_of_watersheds
            ):
                msg = "Length of starting_nodes must equal the" "number_of_watersheds!"
                raise ValueError(msg)
        else:
            large_outlet_ids = grid.boundary_nodes[
                np.argsort(self._stopping_field[grid.boundary_nodes])
            ]
            if number_of_watersheds is None:
                big_enough_watersheds = (
                    self._stopping_field[large_outlet_ids] > outlet_threshold
                )
                starting_nodes = large_outlet_ids[big_enough_watersheds]
            else:
                starting_nodes = large_outlet_ids[-number_of_watersheds:]

        starting_da = self._stopping_field[starting_nodes]
        starting_nodes = np.asarray(starting_nodes)
        if np.any(starting_da < self.threshold) or starting_nodes.size == 0:
            msg = (
                "The number of watersheds requested by the ChannelProfiler is "
                "greater than the number in the domain with sufficent drainage"
                "area."
            )
            raise ValueError(msg)

        self._starting_nodes = starting_nodes

    @property
    def profile_structure(self):
        """
        profile_structure, the channel segment datastructure.
                profile structure is a list of length number_of_watersheds.
                Each element of profile_structure is itself a list of length
                number of stream segments that drain to each of the starting
                nodes. Each stream segment list contains the node ids of a
                stream segment from downstream to upstream.
        """
        return self._profile_structure

    @property
    def distances_upstream(self):
        """
        distances_upstream, the channel segment datastructure.
                A datastructure that parallels profile_structure but holds
                distances upstream instead of node IDs.

                Both lists are number_of_watersheds long.
        """
        return self._distance_along_profile

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
            supplying_nodes = supplying_nodes[np.where(supplying_nodes != i)]

            # if only adding the biggest channel, continue upstream choosing the
            # largest node until no more nodes remain.
            if self._main_channel_only:
                max_drainage = np.argmax(self._stopping_field[supplying_nodes])

                if self._stopping_field[supplying_nodes[max_drainage]] < self.threshold:
                    nodes_to_process = []
                    channel_upstream = False
                else:
                    j = supplying_nodes[max_drainage]

            # if considering multiple channel segments, continue upstream until
            # there are two or more donors with sufficient discharge, then
            # break, returning those nodes as starting points.
            else:

                # get all upstream drainage areas
                upstream_das = self._stopping_field[supplying_nodes]

                # if no nodes upstream exceed the threshold, exit
                if np.sum(upstream_das > self.threshold) == 0:
                    nodes_to_process = []
                    channel_upstream = False

                # otherwise
                else:
                    # if only one upstream node exceeds the threshold, proceed
                    # up the channel.
                    if np.sum(upstream_das > self.threshold) == 1:
                        max_drainage = np.argmax(self._stopping_field[supplying_nodes])
                        j = supplying_nodes[max_drainage]
                    # otherwise provide the multiple upstream nodes to be
                    # processed into a new channel.
                    else:
                        nodes_to_process = supplying_nodes[
                            upstream_das > self.threshold
                        ]
                        channel_upstream = False

        return (channel_segment, nodes_to_process)

    def _create_profile_structure(self):
        """
        Create the profile_IDs data structure for channel network.

        The bound attribute self._profile structure is the channel segment
        datastructure. profile structure is a list of length
        number_of_watersheds. Each element of profile_structure is itself a
        list of length number of stream segments that drain to each of the
        starting nodes. Each stream segment list contains the node ids of a
        stream segment from downstream to upstream.
        """
        self._profile_structure = []

        if self._main_channel_only:
            for i in self._starting_nodes:
                channel_network = []
                (channel_segment, nodes_to_process) = self._get_channel_segment(i)
                channel_network.append(np.array(channel_segment))
                self._profile_structure.append(channel_network)

        else:
            for i in self._starting_nodes:
                channel_network = []
                queue = [i]
                while len(queue) > 0:
                    node_to_process = queue.pop(0)
                    (channel_segment, nodes_to_process) = self._get_channel_segment(
                        node_to_process
                    )
                    channel_network.append(np.array(channel_segment))
                    queue.extend(nodes_to_process)
                self._profile_structure.append(channel_network)

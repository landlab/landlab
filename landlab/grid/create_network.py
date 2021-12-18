# -*- coding: utf-8 -*-
"""
Created on Tue Feb  2 17:50:09 2021

@author: sahrendt
"""
from dataclasses import dataclass
from itertools import tee

import numpy as np
import numpy.typing as npt

from ..components import ChannelProfiler, FlowAccumulator
from ..graph import NetworkGraph
from .network import NetworkModelGrid

try:
    from itertools import pairwise
except ImportError:

    def pairwise(iterable):
        # pairwise('ABCDEFG') --> AB BC CD DE EF FG
        a, b = tee(iterable)
        next(b, None)
        return zip(a, b)


def create_network_from_raster(
    rmg,
    min_channel_thresh=10000,
    outlet_nodes=None,
    method="variable",
    n_widths=20,
    a=9.68,
    b=0.32,
    d_node_spacing=None,
    fields=None,
):

    """Create a NetworkModelGrid from a RasterModelGrid. Default behavior
    spaces nodes a certain number of local channel widths apart across the
    network. If method='constant' is specified, the d_node_spacing value
    is used to space nodes a constant distance apart across the network.

    Parameters
    ----------
    rmg : RasterModelGrid object
        A raster grid used to create a network grid
    min_channel_thresh : float, optional
        Value to use for the minimum drainage area associated with a
        plotted channel segment from the ChannelProfiler. Default values 10000.
    outlet_nodes : Single int value in iterable form, optional
        Iterable containing the node ID of nodes to start the channel profiles
        from in ChannelProfiler. If not provided, the default is the node ID on
        the model grid boundary with the largest terminal drainage area.
    method : string, 'variable' or 'constant'
        Specifies node-spacing method. 'variable' will dynamically update
        node spacing depending on local channel width. 'constant' will use the
        specified 'node_spacing' value to space nodes evenly across the
        network.
    n_widths : float, optional
        Multiplier to define node spacing as a function of local channel width.
    a : float, optional
        Parameter to be used in the Frasson et al. 2019 (GRL) power
        relationship between drainage area and channel width W=a*A^b. Default
        is value is 9.68 from Frassen et. al 2019
    b : float, optional
        Parameter to be used in the Frasson et al. 2019 (GRL) power
        relationship between drainage area and channel width W=a*A^b. Default
        is value is 0.32 from Frassen et. al 2019
    d_node_spacing : float, optional
        Distance value for a constant node spacing along channel segments.
        Must be provided if method is 'constant'.
    fields : iterable, optional
        .at_node fields to map from RasterModelGrid to NetworkModelGrid.
        Formatted as strings inside an iterable object

    Returns
    -------
    NetworkModelGrid object with .at_node['rmg_node_value'] attribute
    listing the RasterModelGrid node ids at each NetworkModelGrid node.

    """

    if "drainage_area" not in rmg.at_node:
        # run flow accumulator for ChannelProfiler
        FlowAccumulator(
            rmg,
            "topographic__elevation",
            flow_director="D8",
            depression_finder="DepressionFinderAndRouter",
        ).run_one_step()

    if "drainage_area" not in rmg.at_node:
        raise ValueError("'drainage_area' field is missing from the grid")

    # delineate channel
    profiler = ChannelProfiler(
        rmg,
        number_of_watersheds=1,
        minimum_channel_threshold=min_channel_thresh,
        outlet_nodes=outlet_nodes,
        main_channel_only=False,
    )
    profiler.run_one_step()

    # obtain watershed key (should only be one)
    wtrshd_key = [k for k in profiler.data_structure.keys()][0]
    # obtain keys for channel segments, keys are raster nodes formatted as
    # tuple for (start, end) for channel segment start and end locations
    channel_segment_keys = profiler.data_structure[wtrshd_key].keys()

    # IDENTIFY CHANNEL SEGMENT CONNECTIVITY -----------------------------------
    # obtain node ids for start and end of every channel segments
    seg_starts = [seg[0] for seg in profiler.data_structure[wtrshd_key].keys()]
    seg_ends = [seg[1] for seg in profiler.data_structure[wtrshd_key].keys()]
    # identify channel connectivity and how to properly link nodes
    # at different channel junctions
    # code does this by identifying the key of the channel seg just downstream
    # and connects first node of upstream channel seg to downstream channel seg
    for seg_key in channel_segment_keys:
        # create empty list to store how segment is connected
        connectivity = []
        connectivity_key = None
        # extract a single segment from the profiler data structure
        seg_i = profiler.data_structure[wtrshd_key][seg_key]
        # ask wehther the start of a segment is the end of another segment
        # (i.e. is this going to be connected downstream?)
        if seg_key[0] in seg_ends:
            connectivity.append("connected downstream")
            # find first segment downstream that should connect to last node of
            # segment upstream
            connect_to_channel_idx = np.argmax(seg_key[0] == seg_ends)
            connectivity_key = (
                seg_starts[connect_to_channel_idx],
                seg_ends[connect_to_channel_idx],
            )
        # ask whether the end of segment is in start of another segment
        # (i.e. is this going to be connected upstream?) if it isn't going to
        # be connected upstream, we will use this to prompt the code to take
        # the last node in the profiler datastructure as the end node. (We
        # don't necessarily need this, if we are okay with some trimming in
        # channel headwaters)
        if seg_key[-1] in seg_starts:
            connectivity.append("connected upstream")
        # note: we do not collect a connectivity_key here, because we
        # connect segments by connecting downstream, not upstream (we don't
        # need to do it twice)

        # lets add this to our segment structure so we can know when and how to
        # connect our segments as we build the grid
        seg_i["connectivity"] = connectivity
        seg_i["connectivity_key"] = connectivity_key

    node_xy = []  # empty list to store paired x,y locations of nmg nodes
    rmg_nodes = (
        []
    )  # empty list to store raster model grid node corresponding to each network model grid node
    links = []  # empty list to store link connections between nodes

    # FUNCTION TO ADD LINKS----------------------------------------------------
    # this is called several times in loop below, hoping it
    # makes testing easier to test once as a function
    def add_link(rmg, all_nodes_xy, all_links, head_node_rmg_id, tail_node_rmg_id):

        """Add link connections to existing list of NetworkModelGrid nodes
        based upon an upstream and downstream RasterModelGrid node id. Also
        checks whether (x, y) values for upstream and downstream nodes exist
        in list of node locations and adds them if necessary.

        Parameters
        ----------
        rmg : RasterModelGrid object
            The RasterModelGrid to which NetworkModelGrid nodes and links will
            be added.
        all_nodes_xy : list of tuples
            List where tuple values for node x and y locations formatted as
            [(x1,y1), (x2,y2)...] already exist or will be stored.
        all_links : list of tuples
            List where tuple values for NetworkModelGrid node ids of upstream
            and downstream nodes for each link already exists or will be
            stored. Formatted as [(id2, id1),(id3, id2)...] where id# '#'
            corresponds to the index of the node entry in all_nodes_xy.
        head_node_rmg_id : int
            Value of the RasterModelGrid node id that corresponds to the
            desired head node of a NetworkModelGrid link.
        tail_node_rmg_id : int
            Value of the RasterModelGrid node id that corresponds to the
            desired head node of a NetworkModelGrid link.

        Returns
        -------
        None.

        """
        # define head node xy value by calling id from raster model grid
        head_node_xy = (
            rmg.x_of_node[head_node_rmg_id],
            rmg.y_of_node[head_node_rmg_id],
        )
        # define a tail node xy value by calling id from raster model grid
        tail_node_xy = (
            rmg.x_of_node[tail_node_rmg_id],
            rmg.y_of_node[tail_node_rmg_id],
        )
        # if these nodes don't already exist in the array of node xy vals from
        # another channel segment, add them
        if head_node_xy not in all_nodes_xy:
            all_nodes_xy.append(head_node_xy)
        if tail_node_xy not in all_nodes_xy:
            all_nodes_xy.append(tail_node_xy)
        # get the index of the head and tail node from our node_xy list
        # this is important to do in case they were added from a previous
        # channel segment. we need to ensure the order of network nodes is
        # correct
        head_node__nmg_id = all_nodes_xy.index(head_node_xy)
        tail_node__nmg_id = all_nodes_xy.index(tail_node_xy)
        # append the head and tail network node ids to the link array
        # this if statement should be sufficient since we only connect links
        # one-way (i.e. downstream) between segments, but if there are still
        # issues with duplicate links something to add would be an additional
        # check for the opposite link combo:
        # if (tail_node__nmg_id, head_node__nmg_id) not in all_links
        if (head_node__nmg_id, tail_node__nmg_id) not in all_links:
            all_links.append((head_node__nmg_id, tail_node__nmg_id))

    # CREATE NETWORK MODEL GRID NODES & LINKS----------------------------------
    # loop over all channel segments and add network model nodes that correspond
    # to a certain cell in the raster model grid
    for seg_key in channel_segment_keys:

        # access data of channel segments
        seg_i = profiler.data_structure[wtrshd_key][seg_key]

        # create list to hold rmg node ids where nmg nodes are located
        nmg_nodes = []

        # identify rmg value of first node in segment
        # first nodes of channel segments will be included in network
        # if they don't already exist
        idx_node = 0
        rmg_node = seg_i["ids"][idx_node]

        # if seg_i is connected dowstream, add link connecting first node of
        # segment to downstream node
        if seg_i["connectivity_key"] is not None:
            channel_key = seg_i["connectivity_key"]
            connecting_seg = profiler.data_structure[wtrshd_key][channel_key]

            # check to make sure there are nmg nodes on downstream segment
            # (there might not be any if the downstream segment is shorter than
            # what the user specifies for the desired network grid node spacing)
            if len(connecting_seg["ids_nmg"]) > 0:
                connect_node = connecting_seg["ids_nmg"][-1]

            # if there are no nmg nodes on the downstream segment
            # it must be too short for calculated node spacing
            # if this is the case, connect upstream segment to first node in
            # dowsntream connecting seg
            else:
                connect_node = connecting_seg["ids"][0]
            # add a link for this connection if necessary
            add_link(
                rmg,
                node_xy,
                links,
                head_node_rmg_id=rmg_node,
                tail_node_rmg_id=connect_node,
            )

        # iterate over segment adding new nodes as long as there are upstream nodes
        # that can be placed on network model grid based upon node spacing
        upstrm_node = True
        while upstrm_node is True:

            # if we haven't already stored the rmg id value for this new node
            # add it to our master list of rmg nodes and sub-list of nmg nodes
            if rmg_node not in rmg_nodes:
                rmg_nodes.append(rmg_node)
                nmg_nodes.append(rmg_node)

            # Assign node spacing as n_channel_widths or a constant value from
            # input params
            if method == "variable":
                # calculate drainage area contributing to this node
                da_node = rmg.at_node["drainage_area"][rmg_node]
                # relate drainage area to river width (convert area to km, width in m)
                # from Frasson et al. 2019 GRL
                w_channel = (a * da_node / (1000 ** 2)) ** b
                # calculate upstream node spacing, n_widths_defines stable node spacing
                node_spacing = n_widths * w_channel
            if method == "constant":
                node_spacing = d_node_spacing

            # if stable node spacing is greater than raster grid resolution
            if node_spacing > rmg.dx:
                # optimal along-channel node location based upon node spacing
                opt_loc = seg_i["distances"][idx_node] + node_spacing
                # define tolerance to not add extra node if opt loc is within half
                # a node spacing away from end of segment
                buffer_tol = 0.5 * node_spacing

                # if we can fit another node on the channel segment
                if opt_loc < (seg_i["distances"][-1] - buffer_tol):
                    # find id of node closest to this location
                    idx_next_node = np.abs(seg_i["distances"] - opt_loc).argmin()
                    # update rmg node with whatever this next node should be
                    rmg_next_node = seg_i["ids"][idx_next_node]

                    # add link from this upstream node to the current node
                    # if necessary
                    add_link(
                        rmg,
                        node_xy,
                        links,
                        head_node_rmg_id=rmg_next_node,
                        tail_node_rmg_id=rmg_node,
                    )

                    # update idx_node and rmg node for next loop
                    rmg_node = rmg_next_node
                    idx_node = idx_next_node

                # if no more nodes can be placed on this segment,
                # move to next segment
                else:
                    upstrm_node = False
                    # add last node in segment to list of node xys
                    last_node_xy = (rmg.x_of_node[rmg_node], rmg.y_of_node[rmg_node])
                    if last_node_xy not in node_xy:
                        node_xy.append(last_node_xy)

            # if no more nodes have stable locations on this segment
            # move to next segment
            else:
                upstrm_node = False

                # check if segment is connected upstream:
                if "connected upstream" in seg_i["connectivity"]:
                    # if we are seeing links on main stem channels that are smaller
                    # then raster model grid resolution, flag this as an error
                    raise ValueError(
                        "main stem link lengths are smaller than grid res."
                        "try increasing n_widths or changing a and b params"
                    )

                # add last node in segment to list of node xys
                last_node_xy = (rmg.x_of_node[rmg_node], rmg.y_of_node[rmg_node])
                if last_node_xy not in node_xy:
                    node_xy.append(last_node_xy)

        # store location of network nodes as raster model ids in channel profiler
        # datastructure. this will be used for joining channel segments later
        seg_i["ids_nmg"] = np.array(nmg_nodes)

    # CREATE NETWORK MODEL GRID OBJECT-----------------------------------------
    x_of_node, y_of_node = zip(*node_xy)

    # Maintain sorting by creating an unsorted network graph and sorting.
    # This process is important to ensure that the fields are assigned to the
    # correct links.
    graph_net = NetworkGraph((y_of_node, x_of_node), links=links, sort=False)
    sorted_nodes, sorted_links, sorted_patches = graph_net.sort()

    # use the sorting information to make a new network model grid.
    nmg = NetworkModelGrid(
        (np.asarray(y_of_node)[sorted_nodes], np.asarray(x_of_node)[sorted_nodes]),
        np.vstack((graph_net.node_at_link_head, graph_net.node_at_link_tail)).T,
    )

    # add RMG node locations and extra fields to network model grid from
    # raster model grid
    nmg.at_node["rmg_node_value"] = np.array(rmg_nodes)[sorted_nodes]
    if fields is None:
        fields = []
    for field in fields:
        nmg.at_node[field] = rmg.at_node[field][nmg.at_node["rmg_node_value"]]

    return nmg


def network_grid_from_raster(
    grid, reducer=None, include="*", exclude=None, minimum_channel_threshold=0.0
):
    """Create a NetworkModelGrid from a RasterModelGrid."""

    if "drainage_area" not in grid.at_node:
        FlowAccumulator(
            grid,
            "topographic__elevation",
            flow_director="D8",
            depression_finder="DepressionFinderAndRouter",
        ).run_one_step()

    if "drainage_area" not in grid.at_node:
        raise ValueError("'drainage_area' field is missing from the grid")

    channel_segments = get_channel_segments(
        grid, minimum_channel_threshold=minimum_channel_threshold
    )

    if reducer is not None:
        channel_segments = [reducer(segment) for segment in channel_segments]

    network_grid = network_grid_from_segments(
        grid, channel_segments, include=include, exclude=exclude
    )

    return network_grid


def get_channel_segments(grid, divergent_okay=False, minimum_channel_threshold=0.0):
    """Extract channel segments from a grid.

    Each segment includes nodes within the segment, upstream segments, and
    downstream segments.

    Parameters
    ----------
    grid : RasterModelGrid
        Grid from which to extract channel segments.
    divergent_okay : bool, optional
        If ``False``, raise an error if the network is divergent.

    Returns
    -------
    segments : list
        Channel segments as lists of grid nodes. Nodes are ordered from downstream
        to upstream.
    """
    # delineate channel
    profiler = ChannelProfiler(
        grid,
        number_of_watersheds=1,
        minimum_channel_threshold=minimum_channel_threshold,
        # outlet_nodes=outlet_nodes,
        main_channel_only=False,
    )
    profiler.run_one_step()

    if len(profiler.data_structure) > 1:
        raise RuntimeError(
            "number of watersheds is greater than the requested number "
            f"({len(profiler.data_structure)} > 1)"
        )

    # obtain watershed key (should only be one)
    watershed = list(profiler.data_structure.keys())[0]

    segments = [
        segment["ids"] for segment in profiler.data_structure[watershed].values()
    ]

    return segments


def network_grid_from_segments(grid, nodes_at_segment, include="*", exclude=None):
    """Create a NetworkModelGrid from channel segments."""
    channel_network = ChannelSegmentConnector(*nodes_at_segment)

    for segment in channel_network.root:
        if len(segment.upstream) == 1:
            if len(segment.upstream[0]) > 0:
                print("segments can be joined")

    xy_of_node = create_xy_of_node(channel_network.root, grid)
    at_node = get_node_fields(
        channel_network.root, grid, include=include, exclude=exclude
    )

    reindex_network_nodes(channel_network.root)
    nodes_at_link = create_network_links(channel_network.root)

    grid = NetworkModelGrid((xy_of_node[:, 1], xy_of_node[:, 0]), links=nodes_at_link)

    for name, values in at_node.items():
        grid.at_node[name] = values[grid._sorted_nodes]

    return grid


class SegmentReducer:
    """Base class for reducing the nodes in a segment."""

    def reduce(self, segment):
        """Reduce the number of nodes in a channel segment."""
        raise NotImplementedError("reduce")

    def __call__(self, segment):
        return self.reduce(segment)


@dataclass
class SpacingAtLeast(SegmentReducer):
    """Remove segment nodes to ensure a minimum along-channel spacing."""

    xy_of_node: npt.ArrayLike
    spacing: npt.ArrayLike = 1.0

    def __post_init__(self):
        self.xy_of_node = np.asarray(self.xy_of_node)
        self.spacing = np.broadcast_to(self.spacing, len(self.xy_of_node))

    def calc_distance_along_segment(self, segment):
        """Calculate the along-channel distance of a segment.

        Parameters
        ----------
        segment : iterable of int
            Indices of nodes along a channel.

        Returns
        -------
        ndarray of int
            Distances to each node along the channel.
        """
        return np.sqrt(
            np.sum(
                np.diff(
                    self.xy_of_node[
                        segment,
                    ],
                    axis=0,
                    prepend=self.xy_of_node[None, segment[0], :],
                )
                ** 2,
                axis=1,
            )
        ).cumsum()

    def reduce(self, segment):
        nodes = _reduce_to_fewest_nodes(
            self.xy_of_node[segment], spacing=self.spacing[segment]
        )
        return np.take(segment, nodes)


@dataclass
class VariableSpacing(SpacingAtLeast):
    """Remove segment nodes to ensure a minimum per-node along-channel spacing."""

    spacing: npt.ArrayLike

    def __post_init__(self):
        self.spacing = np.broadcast_to(self.spacing, len(self.xy_of_node))

    def reduce(self, segment):
        nodes = _reduce_nodes(
            self.calc_distance_along_segment(segment),
            spacing=self.spacing[
                segment,
            ],
        )
        return np.take(segment, nodes)


class JustEndNodes(SegmentReducer):
    """Reduce a segment to just its end nodes."""

    def reduce(self, segment):
        return np.asarray([segment[0], segment[-1]])


@dataclass
class AtMostNodes(SegmentReducer):
    """Reduce a segment to a maximum number of nodes."""

    count: int = 3

    def __post_init__(self):
        if self.count < 2:
            raise ValueError(
                f"unable to reduce a segment to less than two nodes ({self.count})"
            )

    def reduce(self, segment: npt.ArrayLike) -> npt.ArrayLike:
        if self.count < len(segment):
            step = len(segment) // (self.count - 2 + 1)
            reduced_segment = np.append(
                segment[: -step + 1 : step][: self.count - 1], segment[-1]
            )
        else:
            reduced_segment = np.asarray(segment)
        return reduced_segment


def spacing_from_drainage_area(
    drainage_area,
    a=9.68,
    b=0.32,
    n_widths=20.0,
):
    """Calculate channel spacing based on upstream drainage area of each node.

    Parameters
    ----------
    drainage_area : number or ndarray
        Upstream drainage area in km ** 2.

    Returns
    -------
    ndarray
        Node spacing in meters.
    """
    return n_widths * (a * drainage_area / (1000 ** 2)) ** b


def _reduce_nodes(distance_along_segment, spacing=1.0):
    """Reduce the number of nodes in a segment based on a minimum spacing."""
    from bisect import bisect_left

    distance_along_segment = np.asarray(distance_along_segment)
    n_nodes = len(distance_along_segment)
    spacing = np.broadcast_to(spacing, n_nodes)

    nodes = []
    head_node = 0
    while head_node < n_nodes - 1:
        nodes.append(head_node)
        distance_to_tail_node = distance_along_segment[head_node] + spacing[head_node]

        #         tail_node = _find_index_to_nearest(
        #             distance_along_segment[head_node + 1 :], distance_to_tail_node
        #         )
        #
        #         head_node = tail_node + head_node + 1

        tail_node = bisect_left(
            distance_along_segment, distance_to_tail_node, lo=head_node + 1
        )
        head_node = tail_node

    if nodes[-1] < n_nodes - 1:
        nodes.append(n_nodes - 1)

    return nodes


def calc_distance_to_point(point, points):
    """Find the euclidian distance between one point and a set of points."""
    return np.sqrt(np.sum((point - points) ** 2, axis=1))


def _reduce_to_fewest_nodes(xy_of_node, spacing=1.0):
    """Reduce to the fewest number of nodes while maintaining a minimum spacing."""
    xy_of_node = np.asarray(xy_of_node)
    n_nodes = len(xy_of_node)
    spacing = np.broadcast_to(spacing, n_nodes)

    nodes = []
    head_node = 0
    while head_node < n_nodes - 1:
        nodes.append(head_node)
        distance_from_head = calc_distance_to_point(
            xy_of_node[head_node], xy_of_node[head_node + 1 :, :]
        )

        try:
            tail_node = (
                np.where(distance_from_head < spacing[head_node])[0][-1] + head_node + 2
            )
        except IndexError:
            tail_node = head_node + 1
        head_node = tail_node

    if nodes[-1] != n_nodes - 1:
        nodes.append(n_nodes - 1)

    return nodes


def _find_index_to_nearest(array, value):
    """Find index to nearest value in a sorted array of values."""
    array = np.asarray(array)
    if value < array[0]:
        return 0
    elif value >= array[-1]:
        return len(array) - 1
    else:
        ind = np.searchsorted(array, value, side="right")
        return ind - 1 + np.abs(array[ind - 1 : ind + 1] - value).argmin()


class SegmentLinkCollector:
    def __init__(self, links=None):
        if links is None:
            self._links = []
        else:
            self._links = list(links)

    def __call__(self, segment):
        try:
            nodes = [segment.downstream._nodes[-1]]
        except AttributeError:
            nodes = [segment._nodes[0]]
        nodes.extend(segment._nodes[1:])
        for head, tail in pairwise(nodes):
            self._links.append((head, tail))

    @property
    def links(self):
        return self._links


class ChannelSegment:

    """A channel segment.

    Parameters
    ----------
    nodes : iterable of int
        The nodes of the channel, listed from downstream to upstream.
    """

    def __init__(self, nodes):
        self._nodes = None
        self._upstream = []
        self._downstream = None
        self.nodes = nodes

    def __iter__(self):
        yield self
        for upstream in self._upstream:
            for segment in upstream:
                yield segment

    def __len__(self):
        return len(self._nodes)

    def __repr__(self):
        return f"ChannelSegment({self._nodes})"

    def for_each(self, func):
        for segment in self:
            func(segment)

    def iter_downstream(self):
        yield self
        try:
            iter_downstream = self.downstream.iter_downstream
        except AttributeError:
            pass
        else:
            for segment in iter_downstream():
                yield segment

    def count_segments(self, direction="upstream"):
        count = 0
        if direction == "upstream":
            iter = self.__iter__
        elif direction == "downstream":
            iter = self.iter_downstream
        else:
            raise ValueError(f"direction not understood ({direction})")
        for _ in iter():
            count += 1
        return count - 1

    @property
    def downstream(self):
        """The downstream segment."""
        return self._downstream

    @downstream.setter
    def downstream(self, segment):
        self._downstream = segment
        segment.add_upstream(self)

    @property
    def upstream(self):
        """The upstream segments."""
        return tuple(self._upstream)

    def add_upstream(self, segment):
        """Add an upstream segment."""
        self._upstream.append(segment)
        segment._downstream = self

    @property
    def downstream_node(self):
        """The most downstream node of the channel segment."""
        return self._nodes[0]

    @property
    def upstream_node(self):
        """The most upstream node of the channel segment."""
        return self._nodes[-1]

    @property
    def nodes(self):
        """The nodes of the segment, from downstream to upstream."""
        return self._nodes

    @nodes.setter
    def nodes(self, nodes):
        self._nodes = np.array(nodes, copy=True)


class DisconnectedSegmentError(Exception):

    """Raise this exception if a channel segment cannot be connected to a network."""

    pass


class ChannelSegmentConnector:

    """Connect channel segment to form a network."""

    def __init__(self, *args):
        """ChannelSegmentConnector(channel1, channel2, ...)"""
        self._root = None
        self._orphans = []
        for segment in args:
            self.add(segment)

    @property
    def root(self):
        """The root (most downstream) channel of the network."""
        return self._root

    @property
    def orphans(self):
        """Channel segments that are not connected to the main network."""
        return tuple(self._orphans)

    def set_root(self, new_root):
        if self._root is None:
            pass
        elif self._root.downstream_node == new_root.upstream_node:
            new_root.add_upstream(self._root)
        else:
            self._orphans.append(self._root)
        self._root = new_root
        return self._root

    def _add_or_raise(self, new_segment):
        """Try to add a segment to the network, raise an error if disconnected."""
        is_orphan = True
        if (
            self._root is None
            or self._root.downstream_node == new_segment.upstream_node
        ):
            self._root = self.set_root(new_segment)
            is_orphan = False
        else:
            for segment in self._root:
                if new_segment.downstream_node == segment.upstream_node:
                    segment.add_upstream(new_segment)
                    is_orphan = False
                    break

        if is_orphan:
            raise DisconnectedSegmentError()
        return self._root

    def _add_orphans(self):
        """Add orphans to the root."""
        orphans = []
        for orphan in self._orphans:
            try:
                self._root = self._add_or_raise(orphan)
            except DisconnectedSegmentError:
                orphans.append(orphan)
        return orphans

    def add(self, new_segment):
        """Add a new segment to the network."""
        if not isinstance(new_segment, ChannelSegment):
            new_segment = ChannelSegment(new_segment)
        try:
            self._root = self._add_or_raise(new_segment)
        except DisconnectedSegmentError:
            is_orphan = True
        else:
            is_orphan = False

        if not is_orphan:
            self._orphans = self._add_orphans()
        else:
            self._orphans.append(new_segment)

    def __repr__(self):
        return f"ChannelConnector({self._root})"


def create_xy_of_node(network, grid):
    """Create an array of coordinates for each node of a channel network."""
    xy_of_node_collector = SegmentNodeCoordinateCollector(grid)
    network.for_each(xy_of_node_collector)
    return np.asarray(xy_of_node_collector.xy_of_node)


class SegmentNodeCoordinateCollector:
    def __init__(self, grid):
        self._grid = grid
        self._xy_of_node = []

    def __call__(self, segment):
        if segment.downstream is None:
            nodes = segment._nodes
        else:
            nodes = segment._nodes[1:]
        self._xy_of_node.extend(
            zip(self._grid.x_of_node[nodes], self._grid.y_of_node[nodes])
        )

    @property
    def xy_of_node(self):
        return self._xy_of_node


def get_node_fields(network, grid, include="*", exclude=None):
    """Get field values for each node of a channel network.

    Parameters
    ----------
    network : ChannelSegment
        A channel network to extract fields for.
    grid : ModelGrid
        Grid from which to extract fields from.
    include : str or list of str, optional
        Patterns to use for including fields.
    exclude : str or list of str, optional
        Patterns to use for excluding fields.

    Returns
    -------
    at_node : dict
        Dictionary of node fields for each node of a channel network.
    """
    if isinstance(include, str):
        include = [include]

    include = [
        pattern if pattern.startswith("at_") else f"at_node:{pattern}"
        for pattern in include
    ]

    node_fields = set()
    for canonical_name in grid.fields(include=include, exclude=exclude):
        dim, name = canonical_name[len("at_") :].split(":")
        dim == "node" and node_fields.add(name)

    at_node = {}
    for name in node_fields:
        field_value_collector = SegmentFieldCollector(grid.at_node[name])
        network.for_each(field_value_collector)
        at_node[name] = np.asarray(field_value_collector.values)

    return at_node


class SegmentFieldCollector:
    def __init__(self, field):
        self._field = field
        self._values = []

    def __call__(self, segment):
        if segment.downstream is None:
            nodes = segment._nodes
        else:
            nodes = segment._nodes[1:]

        self._values.extend(self._field[nodes])

    @property
    def values(self):
        return self._values


def reindex_network_nodes(network):
    """Reindex the nodes of a channel network."""
    node_reindexer = SegmentNodeReindexer()
    network.for_each(node_reindexer)

    return network


class SegmentNodeReindexer:
    def __init__(self, nodes=None):
        if nodes is None:
            self._nodes = []
        else:
            self._nodes = list(nodes)

    def __call__(self, segment):
        try:
            start = self._nodes[-1] + 1
        except IndexError:
            start = 0

        try:
            downstream_node = segment.downstream._nodes[-1]
        except AttributeError:
            segment._nodes = list(range(start, start + len(segment)))
        else:
            segment._nodes = [downstream_node] + list(
                range(start, start + len(segment) - 1)
            )

        self._nodes.extend(segment._nodes)

    @property
    def nodes(self):
        return self._nodes


def create_network_links(network):
    """Create links between nodes of the channels of a network.

    Parameters
    ----------
    network : ChannelSegment
        Channel network to create links for.

    Returns
    -------
    links : list of tuple
        Links for network as *(head_node, tail_node)*."""
    collect_segment_links = SegmentLinkCollector()
    network.for_each(collect_segment_links)

    return collect_segment_links.links

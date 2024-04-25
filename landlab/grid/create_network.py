"""
Created on Tue Feb  2 17:50:09 2021

@author: sahrendt
"""

from dataclasses import dataclass
from itertools import tee

import numpy as np
import numpy.typing as npt

from ..components import ChannelProfiler
from ..components import FlowAccumulator
from .network import NetworkModelGrid

try:
    from itertools import pairwise
except ImportError:
    # Not available before Python 3.10
    def pairwise(iterable):
        # pairwise('ABCDEFG') --> AB BC CD DE EF FG
        a, b = tee(iterable)
        next(b, None)
        return zip(a, b)


def network_grid_from_raster(
    grid, reducer=None, include="*", exclude=None, minimum_channel_threshold=0.0
):
    """Create a NetworkModelGrid from a RasterModelGrid.

    The optional *reducer* keyword is use to pass a function that reduces the number
    of nodes in each channel segment. The default is to keep all segment nodes.

    Parameters
    ----------
    grid : RasterModelGrid object
        A raster grid used to create a network grid
    reducer : func, optional
        A function used to reduce the number of nodes in each segment. The
        default is to include all segment nodes in the created NetworkModelGrid.
    include : str, or iterable of str, optional
        Glob-style pattern for field names to include.
    exclude : str, or iterable of str, optional
        Glob-style pattern for field names to exclude.
    minimum_channel_threshold : float, optional
        Value to use for the minimum drainage area associated with a
        plotted channel segment from the ChannelProfiler. Default values 10000.

    Returns
    -------
    network : NetworkModelGrid
        NetworkModelGrid object with .at_node['rmg_node_value'] attribute
        listing the RasterModelGrid node ids at each NetworkModelGrid node.
    """

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
        If ``False``, raise an error if the network is divergent (i.e. a channel
        segment has more than one downstream segments).

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
        if len(segment.upstream) == 1 and len(segment.upstream[0]) > 0:
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
                    self.xy_of_node[segment,],
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
class AlongChannelSpacingAtLeast(SpacingAtLeast):
    """Remove segment nodes to ensure a minimum per-node along-channel spacing."""

    def reduce(self, segment):
        nodes = _reduce_nodes(
            self.calc_distance_along_segment(segment),
            spacing=self.spacing[segment,],
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
    return n_widths * (a * drainage_area / (1000**2)) ** b


def _reduce_nodes(distance_along_segment, spacing=1.0):
    """Reduce the number of nodes in a segment based on a minimum spacing.

    Parameters
    ----------
    distance_along_segment : array of float
        Distance along a segment to each of the segment's nodes.
    spacing : float or array of float, optional
        Minimum spacing for each node along a segment. If a scalar,
        a constant spacing is used along the segment.

    Returns
    -------
    list : nodes
        Indices of nodes to retain after the reduction.

    Examples
    --------
    >>> from landlab.grid.create_network import _reduce_nodes

    Maintain a spacing of at least 1.75.

    >>> distance = [0.0, 1.0, 2.0, 3.0, 4.0]
    >>> _reduce_nodes(distance, spacing=1.75)
    [0, 2, 4]

    If the requested minimum spacing is already smaller than the
    initial spacing, keep all the nodes.

    >>> distance = [0.0, 1.0, 2.0, 3.0, 4.0]
    >>> _reduce_nodes(distance, spacing=0.5)
    [0, 1, 2, 3, 4]

    The spacing can be variable from node to node.

    >>> distance = [0.0, 1.0, 2.0, 3.0, 4.0]
    >>> _reduce_nodes(distance, spacing=[0.5, 1.0, 2.0, 1.0, 0.5])
    [0, 1, 2, 4]

    The end nodes are always retained.

    >>> distance = [0.0, 1.0, 2.0, 3.0, 4.0]
    >>> _reduce_nodes(distance, spacing=100.0)
    [0, 4]

    """
    from bisect import bisect_left

    distance_along_segment = np.asarray(distance_along_segment)
    n_nodes = len(distance_along_segment)
    spacing = np.broadcast_to(spacing, n_nodes)

    nodes = []
    head_node = 0
    while head_node < n_nodes - 1:
        nodes.append(head_node)
        distance_to_tail_node = distance_along_segment[head_node] + spacing[head_node]

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
    """Reduce to the fewest number of nodes while maintaining a minimum spacing.

    Parameters
    ----------
    xy_of_node : array of float shape (n_nodes, 2)
        x and y coordinates of each node along a segment.
    spacing : float or array of float, optional
        Minimum spacing for each node along a segment. If a scalar,
        a constant spacing is used along the segment.

    Returns
    -------
    list : nodes
        Indices of nodes to retain after the reduction.

    >>> from landlab.grid.create_network import _reduce_to_fewest_nodes

    Maintain a spacing of at least 1.75.

    >>> xy_of_node = [[0.0, 0.0], [1.0, 0.0], [2.0, 0.0], [3.0, 0.0], [4.0, 0.0]]
    >>> _reduce_to_fewest_nodes(xy_of_node, spacing=1.75)
    [0, 2, 4]

    If the requested minimum spacing is already smaller than the
    initial spacing, keep all the nodes.

    >>> xy_of_node = [[0.0, 0.0], [1.0, 0.0], [2.0, 0.0], [3.0, 0.0], [4.0, 0.0]]
    >>> _reduce_to_fewest_nodes(xy_of_node, spacing=0.5)
    [0, 1, 2, 3, 4]

    The spacing can be variable from node to node.

    >>> xy_of_node = [[0.0, 0.0], [1.0, 0.0], [2.0, 0.0], [3.0, 0.0], [4.0, 0.0]]
    >>> _reduce_to_fewest_nodes(xy_of_node, spacing=[0.5, 1.0, 2.0, 1.0, 0.5])
    [0, 1, 2, 4]

    The end nodes are always retained.

    >>> xy_of_node = [[0.0, 0.0], [1.0, 0.0], [2.0, 0.0], [3.0, 0.0], [4.0, 0.0]]
    >>> _reduce_to_fewest_nodes(xy_of_node, spacing=100.0)
    [0, 4]
    """
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


class SegmentLinkCollector:
    """Collect links between nodes of segments."""

    def __init__(self, links=None):
        if links is None:
            self._links = []
        else:
            self._links = list(links)

    def __call__(self, segment):
        """Add links between segment nodes to those previously collected."""
        try:
            nodes = [segment.downstream._nodes[-1]]
        except AttributeError:
            nodes = [segment._nodes[0]]
        nodes.extend(segment._nodes[1:])
        for head, tail in pairwise(nodes):
            self._links.append((head, tail))

    @property
    def links(self):
        """Head and tail nodes of all collected links."""
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
            yield from upstream

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
            yield from iter_downstream()

    def count_segments(self, direction="upstream"):
        # count = 0
        if direction == "upstream":
            iter = self.__iter__
        elif direction == "downstream":
            iter = self.iter_downstream
        else:
            raise ValueError(f"direction not understood ({direction})")
        # for _ in iter():
        #     count += 1
        # return count - 1
        return sum(1 for _ in iter()) - 1

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
    """Connect channel segments to form a network."""

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
    """Collect xy coordinates for each node along segments."""

    def __init__(self, grid):
        self._grid = grid
        self._xy_of_node = []

    def __call__(self, segment):
        """Add coordinates of the nodes of a segment to previously collected coordinates."""
        if segment.downstream is None:
            nodes = segment._nodes
        else:
            nodes = segment._nodes[1:]
        self._xy_of_node.extend(
            zip(self._grid.x_of_node[nodes], self._grid.y_of_node[nodes])
        )

    @property
    def xy_of_node(self):
        """Coordinates of all collected nodes."""
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
    """Collect field values for each node along segments."""

    def __init__(self, field):
        self._field = field
        self._values = []

    def __call__(self, segment):
        """Add field values for nodes along a segment to previously collected values."""
        if segment.downstream is None:
            nodes = segment._nodes
        else:
            nodes = segment._nodes[1:]

        self._values.extend(self._field[nodes])

    @property
    def values(self):
        """Field values of all collected nodes."""
        return self._values


def reindex_network_nodes(network):
    """Reindex the nodes of a channel network."""
    node_reindexer = SegmentNodeReindexer()
    network.for_each(node_reindexer)

    return network


class SegmentNodeReindexer:
    """Reindex nodes along segments."""

    def __init__(self, nodes=None):
        if nodes is None:
            self._nodes = []
        else:
            self._nodes = list(nodes)

    def __call__(self, segment):
        """Reindex nodes of a segment based on previously collected nodes."""
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
        """Reindexed nodes of all collected nodes."""
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

import hypothesis.extra.numpy as hynp
import numpy as np
import pytest
from hypothesis import given
from hypothesis.strategies import composite
from hypothesis.strategies import floats
from hypothesis.strategies import integers
from hypothesis.strategies import lists
from numpy.testing import assert_array_equal

from landlab import RasterModelGrid
from landlab.grid.create_network import AlongChannelSpacingAtLeast
from landlab.grid.create_network import AtMostNodes
from landlab.grid.create_network import ChannelSegment
from landlab.grid.create_network import ChannelSegmentConnector
from landlab.grid.create_network import JustEndNodes
from landlab.grid.create_network import SegmentLinkCollector
from landlab.grid.create_network import SegmentNodeCoordinateCollector
from landlab.grid.create_network import SegmentNodeReindexer
from landlab.grid.create_network import SpacingAtLeast
from landlab.grid.create_network import _reduce_nodes
from landlab.grid.create_network import _reduce_to_fewest_nodes
from landlab.grid.create_network import create_network_links
from landlab.grid.create_network import create_xy_of_node
from landlab.grid.create_network import get_node_fields
from landlab.grid.create_network import network_grid_from_raster
from landlab.grid.create_network import network_grid_from_segments
from landlab.grid.create_network import pairwise
from landlab.grid.create_network import reindex_network_nodes
from landlab.grid.create_network import spacing_from_drainage_area


@given(
    drainage_area=hynp.arrays(
        dtype=hynp.floating_dtypes(),
        shape=hynp.array_shapes(),
        elements=floats(min_value=0, width=16),
    )
)
def test_calc_spacing_always_positive(drainage_area):
    assert np.all(spacing_from_drainage_area(drainage_area) >= 0.0)


@given(
    drainage_area=hynp.arrays(
        dtype=hynp.floating_dtypes(),
        shape=hynp.array_shapes(),
        elements=floats(min_value=0, width=16),
    )
)
def test_calc_spacing_unit_keywords(drainage_area):
    spacing = spacing_from_drainage_area(drainage_area, a=1, b=1, n_widths=1)

    assert np.allclose(spacing, drainage_area / 1e6)


@given(nodes=lists(integers(), min_size=2, max_size=1024))
def test_channel_segment(nodes):
    segment = ChannelSegment(nodes)
    assert segment.downstream_node == nodes[0]
    assert segment.upstream_node == nodes[-1]
    assert_array_equal(segment.nodes, nodes)
    assert len(segment) == len(nodes)


@given(nodes=lists(integers(), min_size=2, max_size=1024))
def test_channel_segment_set_nodes(nodes):
    segment = ChannelSegment([0, 1])

    segment.nodes = nodes
    assert segment.downstream_node == nodes[0]
    assert segment.upstream_node == nodes[-1]
    assert_array_equal(segment.nodes, nodes)
    assert len(segment) == len(nodes)


def test_channel_segment_add_downstream_node():
    segment = ChannelSegment([0, 1])
    downstream = ChannelSegment([5, 6])
    assert segment.downstream is None
    assert len(downstream.upstream) == 0

    segment.downstream = downstream

    assert segment.downstream is downstream
    assert segment in downstream.upstream


def test_channel_segment_add_upstream_node():
    segment = ChannelSegment([0, 1])
    upstream = ChannelSegment([5, 6])
    assert len(segment.upstream) == 0
    assert upstream.downstream is None

    segment.add_upstream(upstream)

    assert upstream in segment.upstream
    assert upstream.downstream is segment


@given(segments=lists(lists(integers(), min_size=2, max_size=1024), min_size=1))
def test_channel_segment_many_upstream(segments):
    segments = [ChannelSegment(segment) for segment in segments]
    root = segments[0]
    for current, next in pairwise(segments):
        current.add_upstream(next)

    assert root.count_segments(direction="upstream") == len(segments) - 1
    assert root.count_segments(direction="downstream") == 0


@given(segments=lists(lists(integers(), min_size=2, max_size=1024), min_size=1))
def test_channel_segment_many_flat_upstream(segments):
    segments = [ChannelSegment(segment) for segment in segments]
    root = segments[0]
    for segment in segments[1:]:
        root.add_upstream(segment)
    assert root.downstream is None
    assert len(root.upstream) == len(segments) - 1

    assert root.count_segments(direction="upstream") == len(segments) - 1
    assert root.count_segments(direction="downstream") == 0


@given(segments=lists(lists(integers(), min_size=2, max_size=1024), min_size=1))
def test_channel_segment_many_downstream(segments):
    segments = [ChannelSegment(segment) for segment in segments]
    root = segments[0]
    for current, next in pairwise(segments):
        current.downstream = next
    root = segments[0]
    leaf = segments[-1]

    assert root.count_segments(direction="upstream") == 0
    assert root.count_segments(direction="downstream") == len(segments) - 1

    assert leaf.count_segments(direction="upstream") == len(segments) - 1
    assert leaf.count_segments(direction="downstream") == 0


@given(nodes=lists(integers(), min_size=2, max_size=1024))
def test_channel_segment_for_each(nodes):
    all_nodes = []

    def collect_nodes(segment):
        all_nodes.extend(list(segment.nodes))

    segment = ChannelSegment(nodes)
    segment.for_each(collect_nodes)

    assert_array_equal(all_nodes, segment.nodes)


def test_connector_add_upstream():
    segment = ChannelSegment([0, 1])
    connector = ChannelSegmentConnector(segment)

    assert connector.root is segment
    assert len(connector.orphans) == 0

    connector.add(ChannelSegment([1, 2]))
    assert connector.root is segment
    assert connector.root.count_segments(direction="upstream") == 1
    assert connector.root.downstream is None


def test_connector_add_downstream():
    segment_1 = ChannelSegment([0, 1])
    segment_2 = ChannelSegment([2, 0])
    connector = ChannelSegmentConnector(segment_1)

    connector.add(segment_2)
    assert connector.root is segment_2
    assert connector.root.count_segments(direction="upstream") == 1
    assert connector.root.downstream is None


def test_connector_add_orphan():
    segment_1 = ChannelSegment([0, 1])
    segment_2 = ChannelSegment([2, 3])
    connector = ChannelSegmentConnector(segment_1)

    connector.add(segment_2)

    assert connector.root is segment_1
    assert connector.root.count_segments(direction="upstream") == 0
    assert connector.root.downstream is None

    assert len(connector.orphans) == 1
    assert connector.orphans == (segment_2,)

    connector.add(ChannelSegment([1, 2]))
    assert connector.root.count_segments(direction="upstream") == 2
    assert connector.orphans == ()


_grid_dims_to_test = integers(min_value=3, max_value=128)


@composite
def shape_and_indices(draw, elements=_grid_dims_to_test):
    shape = draw(lists(elements, min_size=2, max_size=2))
    indices = draw(
        lists(
            integers(min_value=0, max_value=shape[0] * shape[1] - 1),
            min_size=1,
            max_size=1024,
        )
    )
    return shape, indices


@given(shape_and_segment=shape_and_indices())
def test_construct_xy_of_node(shape_and_segment):
    shape, segment = shape_and_segment
    grid = RasterModelGrid(shape)

    collect_coordinates = SegmentNodeCoordinateCollector(grid)
    collect_coordinates(ChannelSegment(segment))
    xy_of_node = collect_coordinates.xy_of_node

    assert len(xy_of_node) == len(segment)

    x_of_node, y_of_node = zip(*xy_of_node)
    assert (x_of_node == grid.x_of_node[segment]).all()
    assert (y_of_node == grid.y_of_node[segment]).all()


@given(nodes=lists(integers(), min_size=0, max_size=1024))
def test_reindex_segment_nodes_orphan(nodes):
    segment = ChannelSegment(nodes)
    reindex = SegmentNodeReindexer()
    reindex(segment)
    assert segment.nodes == list(range(len(nodes)))


@given(nodes=lists(integers(), min_size=0, max_size=1024), last_node=integers())
def test_reindex_segment_nodes_with_last_node(nodes, last_node):
    segment = ChannelSegment(nodes)
    reindex = SegmentNodeReindexer(nodes=[last_node])
    reindex(segment)
    assert segment.nodes == list(range(last_node + 1, last_node + 1 + len(nodes)))


@given(nodes=lists(integers(), min_size=0, max_size=1024), last_node=integers())
def test_reindex_segment_nodes_with_downstream(nodes, last_node):
    root = ChannelSegment([0, 1])
    segment = ChannelSegment(nodes)
    root.add_upstream(segment)

    reindex = SegmentNodeReindexer(nodes=[last_node])
    reindex(segment)

    assert segment.nodes[0] == root.nodes[-1]
    assert segment.nodes[1:] == list(range(last_node + 1, last_node + len(nodes)))


@given(nodes=lists(integers(), min_size=2, max_size=1024))
def test_create_links(nodes):
    segment = ChannelSegment(nodes)

    collect_links = SegmentLinkCollector()
    collect_links(segment)
    links = collect_links.links

    assert len(links) == len(segment) - 1
    heads, tails = zip(*links)
    assert list(heads) == nodes[:-1]
    assert list(tails) == nodes[1:]


@given(nodes=lists(integers(), min_size=2, max_size=1024))
def test_create_links_with_existing(nodes):
    segment = ChannelSegment(nodes)

    collect_links = SegmentLinkCollector(links=[(1, 2), (3, 4)])
    collect_links(segment)
    links = collect_links.links

    assert links[:2] == [(1, 2), (3, 4)]
    assert len(links[2:]) == len(segment) - 1
    heads, tails = zip(*links[2:])
    assert list(heads) == nodes[:-1]
    assert list(tails) == nodes[1:]


@given(nodes=lists(integers(), min_size=2, max_size=1024))
def test_create_links_with_downstream(nodes):
    root = ChannelSegment([0, 1])
    segment = ChannelSegment(nodes)
    root.add_upstream(segment)

    collect_links = SegmentLinkCollector()
    collect_links(segment)
    links = collect_links.links

    assert len(links) == len(segment) - 1
    assert links[0] == (root.nodes[-1], segment.nodes[1])

    if len(links) > 1:
        heads, tails = zip(*links[1:])
        assert list(heads) == nodes[1:-1]
        assert list(tails) == nodes[2:]


def test_reindex_network_nodes():
    root = ChannelSegmentConnector([10, 11, 12], [12, 13], [12, 14], [14, 15]).root
    reindex_network_nodes(root)
    assert list(root.nodes) == [0, 1, 2]
    assert list(root.upstream[0].nodes) == [2, 3]
    assert list(root.upstream[1].nodes) == [2, 4]
    assert list(root.upstream[1].upstream[0].nodes) == [4, 5]


def test_create_network_links():
    root = ChannelSegmentConnector([0, 1, 2], [2, 3], [2, 4], [4, 5]).root
    links = create_network_links(root)
    assert links == [(0, 1), (1, 2), (2, 3), (2, 4), (4, 5)]


def test_graph_from_segments():
    r"""
    ::
      *
      |
      *   *
       \ /
    *   *
    |   |
    *   *     *
    |   |    /
    *   * - *
     \ /
      *
       \
        *
        |
        *
    """
    grid = RasterModelGrid((8, 6))
    grid.at_node["z"] = list(range(grid.number_of_nodes))
    segments = [
        [3, 9, 14],
        [14, 19, 25, 31],
        [14, 21],
        [21, 27, 33],
        [33, 40],
        [33, 38, 44],
        [21, 22, 29],
    ]
    graph = network_grid_from_segments(grid, segments)
    assert graph.number_of_nodes == 14
    assert graph.number_of_links == 13

    assert list(zip(graph.x_of_node, graph.y_of_node)) == [
        (3.0, 0.0),
        (3.0, 1.0),
        (2.0, 2.0),
        (1.0, 3.0),
        (3.0, 3.0),
        (4.0, 3.0),
        (1.0, 4.0),
        (3.0, 4.0),
        (5.0, 4.0),
        (1.0, 5.0),
        (3.0, 5.0),
        (2.0, 6.0),
        (4.0, 6.0),
        (2.0, 7.0),
    ]

    assert "z" in graph.at_node
    assert list(graph.at_node["z"]) == [
        3,
        9,
        14,
        19,
        21,
        22,
        25,
        27,
        29,
        31,
        33,
        38,
        40,
        44,
    ]


def test_reduce_nodes():
    nodes = _reduce_nodes([0.0, 1.0, 2.0, 3.0, 3.5], spacing=1.0)
    assert nodes == [0, 1, 2, 3, 4]

    nodes = _reduce_nodes([0.0, 1.0, 2.0, 3.0, 3.5], spacing=0.5)
    assert nodes == [0, 1, 2, 3, 4]

    nodes = _reduce_nodes([0.0, 1.0, 2.0, 3.0, 4.0], spacing=1.75)
    assert nodes == [0, 2, 4]

    nodes = _reduce_nodes([0.0, 1.0, 2.0, 3.0, 4.0], spacing=[1.0, 1.0, 2.0, 2.0, 2.0])
    assert nodes == [0, 1, 2, 4]

    nodes = _reduce_nodes([0.0, 1.0, 2.0, 3.0, 4.0], spacing=1000.0)
    assert nodes == [0, 4]

    nodes = _reduce_nodes([0.0, 1.0, 2.0, 3.0, 4.0, 5.0], spacing=2.0)
    assert nodes == [0, 2, 4, 5]


def test_reduce_to_fewest_nodes():
    x = [0.0, 1.0, 2.0, 3.0, 3.5]
    y = [0.0] * len(x)
    nodes = _reduce_to_fewest_nodes(list(zip(x, y)), spacing=1.0)
    assert nodes == [0, 1, 2, 3, 4]

    nodes = _reduce_to_fewest_nodes(list(zip(x, y)), spacing=0.5)
    assert nodes == [0, 1, 2, 3, 4]

    x = [0.0, 1.0, 2.0, 3.0, 4.0]
    y = [0.0] * len(x)
    nodes = _reduce_to_fewest_nodes(list(zip(x, y)), spacing=1.75)
    assert nodes == [0, 2, 4]

    nodes = _reduce_to_fewest_nodes(list(zip(x, y)), spacing=[1.0, 1.0, 2.0, 2.0, 2.0])
    assert nodes == [0, 1, 2, 4]

    nodes = _reduce_to_fewest_nodes(list(zip(x, y)), spacing=1000.0)
    assert nodes == [0, 4]

    x = [0.0, 1.0, 2.0, 3.0, 4.0, 5.0]
    y = [0.0] * len(x)
    nodes = _reduce_to_fewest_nodes(list(zip(x, y)), spacing=2.0)
    assert nodes == [0, 2, 4, 5]


def test_reduce_nodes_stay_the_same():
    nodes = _reduce_nodes([0.0, 1.0, 2.0, 3.0, 4.0], spacing=1.0)
    assert nodes == [0, 1, 2, 3, 4]

    nodes = _reduce_nodes([0.0, 2.0, 4.0, 6.0, 8.0], spacing=2.0)
    assert nodes == [0, 1, 2, 3, 4]

    nodes = _reduce_nodes([0.0, 1.0, 2.0, 3.0, 4.0], spacing=0.5)
    assert nodes == [0, 1, 2, 3, 4]

    nodes = _reduce_nodes([0.0, 1.0, 3.0, 6.0, 10.0], spacing=[1, 2, 3, 4, 5])
    assert nodes == [0, 1, 2, 3, 4]


@pytest.mark.parametrize(
    "x,spacing",
    [
        ([0.0, 1.0, 2.0, 3.0, 4.0], 1.0),
        ([0.0, 2.0, 4.0, 6.0, 8.0], 2.0),
        ([0.0, 1.0, 2.0, 3.0, 4.0], 0.5),
        ([0.0, 1.0, 3.0, 6.0, 10.0], [1, 2, 3, 4, 5]),
    ],
)
def test_reduce_to_fewest_nodes_stay_the_same(x, spacing):
    y = [0.0] * len(x)
    nodes = _reduce_to_fewest_nodes(list(zip(x, y)), spacing=spacing)
    assert nodes == [0, 1, 2, 3, 4]


@given(
    spacing=hynp.arrays(
        dtype=float,
        shape=hynp.array_shapes(min_dims=1, max_dims=1, min_side=2),
        elements=floats(min_value=1e-3, max_value=1e3),
    )
)
def test_reduce_nodes_min_max_spacing(spacing):
    distance_along_segment = np.cumsum(spacing)

    if np.any(np.diff(distance_along_segment) <= 0):
        raise ValueError(f"array not sorted ({distance_along_segment})")

    nodes = _reduce_nodes(distance_along_segment, spacing=spacing.min())
    assert np.all(nodes == np.arange(len(spacing)))

    nodes = _reduce_nodes(
        distance_along_segment,
        spacing=distance_along_segment[-1] - distance_along_segment[0],
    )
    assert nodes == [0, len(spacing) - 1]


@given(
    spacing=hynp.arrays(
        dtype=float,
        shape=hynp.array_shapes(min_dims=1, max_dims=1, min_side=2),
        elements=floats(min_value=1e-3, max_value=1e3),
    )
)
def test_reduce_to_fewest_nodes_min_max_spacing(spacing):
    distance_along_segment = np.cumsum(spacing)

    if np.any(np.diff(distance_along_segment) <= 0):
        raise ValueError(f"array not sorted ({distance_along_segment})")

    xy_of_node = list(zip(distance_along_segment, [0.0] * len(distance_along_segment)))
    min_spacing = np.diff(distance_along_segment).min()

    nodes = _reduce_to_fewest_nodes(xy_of_node, spacing=min_spacing)
    assert np.all(nodes == np.arange(len(spacing)))

    nodes = _reduce_to_fewest_nodes(
        xy_of_node,
        spacing=distance_along_segment[-1] - distance_along_segment[0],
    )
    assert nodes == [0, len(spacing) - 1]


def test_educe_to_fewest_nodes_wraparound():
    x = [0.0, 0.0, 0.0, 1.0, 1.0, 1.0]
    y = [0.0, 1.0, 2.0, 2.0, 1.0, 0.0]

    assert _reduce_to_fewest_nodes(list(zip(x, y)), spacing=1.001) == [0, 5]


def test_create_xy_of_node_with_branch():
    grid = RasterModelGrid((3, 4), xy_spacing=(2.0, 3.0))
    network = ChannelSegmentConnector([4, 5], [5, 2, 3], [5, 10, 11])
    xy_of_node = create_xy_of_node(network.root, grid)

    assert np.allclose(xy_of_node, [[0, 3], [2, 3], [4, 0], [6, 0], [4, 6], [6, 6]])


@given(
    nodes=lists(
        integers(min_value=0, max_value=1023),
        min_size=2,
        max_size=1024,
    ),
)
def test_create_xy_of_node_one_segement(nodes):
    grid = RasterModelGrid((16, 64), xy_spacing=(2.0, 3.0))
    network = ChannelSegmentConnector(nodes)
    xy_of_node = create_xy_of_node(network.root, grid)

    assert np.allclose(xy_of_node[:, 0], grid.x_of_node[nodes])
    assert np.allclose(xy_of_node[:, 1], grid.y_of_node[nodes])


def test_xy_of_node_if_not_network_root():
    grid = RasterModelGrid((3, 4), xy_spacing=(2.0, 3.0))
    network = ChannelSegmentConnector([4, 5], [5, 2], [5, 10, 11], [2, 3], [2, 7])

    base = network.root.upstream[0]
    xy_of_node = create_xy_of_node(base, grid)

    assert np.allclose(xy_of_node, [[4, 0], [6, 0], [6, 3]])


def test_get_node_fields_one_field():
    grid = RasterModelGrid((3, 4))
    grid.at_node["foo"] = np.arange(12) * 10
    network = ChannelSegmentConnector([0, 5], [5, 6, 7], [5, 9])

    fields = get_node_fields(network.root, grid)
    assert list(fields) == ["foo"]
    assert_array_equal(fields["foo"], [0, 50, 60, 70, 90])


def test_get_node_fields_two_fields():
    grid = RasterModelGrid((3, 4))
    grid.at_node["foo"] = np.arange(12) * 10
    grid.at_node["bar"] = np.arange(12) * 100

    network = ChannelSegmentConnector([0, 5], [5, 6, 7], [5, 9])

    fields = get_node_fields(network.root, grid)
    assert sorted(fields) == ["bar", "foo"]
    assert_array_equal(fields["foo"], [0, 50, 60, 70, 90])
    assert_array_equal(fields["bar"], [0, 500, 600, 700, 900])


def test_get_node_fields_include():
    grid = RasterModelGrid((3, 4))
    grid.at_node["foo"] = np.arange(12) * 10
    grid.at_node["bar"] = np.arange(12) * 100
    grid.at_node["baz"] = np.arange(12) * 1000

    network = ChannelSegmentConnector([0, 5], [5, 6, 7], [5, 9])

    fields = get_node_fields(network.root, grid, include="at_node:f*")
    assert list(fields) == ["foo"]
    assert_array_equal(fields["foo"], [0, 50, 60, 70, 90])

    fields = get_node_fields(network.root, grid, include="at_node:b*")
    assert sorted(fields) == ["bar", "baz"]
    assert_array_equal(fields["bar"], [0, 500, 600, 700, 900])
    assert_array_equal(fields["baz"], [0, 5000, 6000, 7000, 9000])


def test_get_node_fields_exclude():
    grid = RasterModelGrid((3, 4))
    grid.add_empty("foo", at="node")
    grid.add_empty("bar", at="node")
    grid.add_empty("baz", at="node")

    network = ChannelSegmentConnector([0, 5], [5, 6, 7], [5, 9])

    expected = get_node_fields(network.root, grid, include="at_node:b*")
    actual = get_node_fields(network.root, grid, exclude="at_node:f*")

    assert actual.keys() == expected.keys()
    for name in actual:
        assert_array_equal(actual[name], expected[name])


def test_get_node_fields_ignore_non_node_fields():
    grid = RasterModelGrid((3, 4))
    grid.add_empty("foo", at="node")
    grid.add_empty("bar", at="node")
    grid.add_empty("baz", at="link")

    network = ChannelSegmentConnector([0, 5], [5, 6, 7], [5, 9])

    fields = get_node_fields(network.root, grid, include="*")

    assert sorted(fields) == ["bar", "foo"]


def test_network_grid_from_raster():
    grid = RasterModelGrid((4, 5))
    grid.at_node["topographic__elevation"] = np.flipud(
        np.asarray(
            [
                [4, 4, 4, 4, 4],
                [4, 4, 2, 4, 4],
                [4, 4, 1, 4, 4],
                [4, 4, 0, 4, 4],
            ],
            dtype=float,
        )
    )

    network = network_grid_from_raster(grid)
    assert network.number_of_nodes == 7
    assert network.number_of_links == 6
    assert np.allclose(
        network.xy_of_node,
        [
            [2.0, 0.0],
            [1.0, 1.0],
            [2.0, 1.0],
            [3.0, 1.0],
            [1.0, 2.0],
            [2.0, 2.0],
            [3.0, 2.0],
        ],
    )
    assert_array_equal(
        network.nodes_at_link, [[0, 2], [1, 2], [2, 3], [4, 2], [2, 5], [2, 6]]
    )


@given(nodes=lists(integers(), min_size=2))
def test_reducer_just_end_nodes(nodes):
    reduce = JustEndNodes()
    assert_array_equal(reduce(nodes), [nodes[0], nodes[-1]])


@given(nodes=lists(integers(), min_size=3))
def test_reducer_min_three_nodes(nodes):
    reduce = AtMostNodes()
    assert_array_equal(reduce(nodes), [nodes[0], nodes[len(nodes) // 2], nodes[-1]])


@given(nodes=lists(integers(), min_size=2))
def test_reducer_min_ndoes_matches_just_end_nodes(nodes):
    just_end_nodes = JustEndNodes()
    at_most_two_nodes = AtMostNodes(count=2)
    assert_array_equal(just_end_nodes(nodes), at_most_two_nodes(nodes))


@given(nodes=lists(integers(), min_size=2))
def test_reducer_min_nodes_no_change(nodes):
    reduce = AtMostNodes(count=len(nodes) + 1)
    assert_array_equal(reduce(nodes), nodes)


@pytest.mark.parametrize("count", [-1, 0, 1])
def test_reducer_min_nodes_less_than_two(count):
    with pytest.raises(ValueError):
        AtMostNodes(count=count)


def test_reducer_spacing_at_least():
    grid = RasterModelGrid((3, 6), xy_spacing=(3.0, 4.0))
    reduce = SpacingAtLeast(xy_of_node=grid.xy_of_node)

    assert_array_equal(
        reduce.calc_distance_along_segment([6, 7, 8, 9, 10, 11]), [0, 3, 6, 9, 12, 15]
    )
    assert_array_equal(reduce.calc_distance_along_segment([0, 7, 14]), [0, 5, 10])
    assert_array_equal(reduce.calc_distance_along_segment([0, 6, 12]), [0, 4, 8])

    reduce = SpacingAtLeast(xy_of_node=grid.xy_of_node, spacing=3.0)
    assert_array_equal(reduce([6, 7, 8, 9, 10, 11]), [6, 7, 8, 9, 10, 11])

    reduce = SpacingAtLeast(xy_of_node=grid.xy_of_node, spacing=1.5)
    assert_array_equal(reduce([6, 7, 8, 9, 10, 11]), [6, 7, 8, 9, 10, 11])

    reduce = SpacingAtLeast(xy_of_node=grid.xy_of_node, spacing=6.0)

    assert_array_equal(reduce([6, 7, 8, 9, 10, 11]), [6, 8, 10, 11])


def test_reducer_spacing_at_least_variable():
    xy_of_node = [[0, 0], [1, 0], [2, 0], [3, 0], [4, 0], [5, 0]]
    spacing = [1, 2, 3, 4, 5, 6]

    reduce = SpacingAtLeast(xy_of_node=xy_of_node, spacing=spacing)

    assert_array_equal(reduce([0, 1, 2, 3, 4, 5]), [0, 1, 3, 5])
    assert_array_equal(reduce([0, 1, 2, 3, 4, 5]), reduce(reduce([0, 1, 2, 3, 4, 5])))


@given(
    xy_of_node=hynp.arrays(
        dtype=float,
        shape=hynp.array_shapes(min_dims=1, max_dims=1, min_side=4, max_side=128),
        elements=integers(min_value=-1024, max_value=1024),
        unique=True,
    ),
    spacing=floats(min_value=0.0, exclude_min=True),
)
def test_reducer_spacing_at_least_all_greater(xy_of_node, spacing):
    xy_of_node = xy_of_node[: len(xy_of_node) - len(xy_of_node) % 2].reshape((-1, 2))
    xy_of_node /= 100.0

    segment = np.arange(len(xy_of_node))
    reduce = SpacingAtLeast(xy_of_node=xy_of_node, spacing=spacing)
    reduced_segment = reduce(segment)

    distance_along_segment = reduce.calc_distance_along_segment(reduced_segment[:-1])

    assert reduced_segment[0] == segment[0]
    assert reduced_segment[-1] == segment[-1]
    assert np.all(np.diff(distance_along_segment) >= spacing)


def test_reducer_along_channel_spacing_at_least_variable():
    xy_of_node = [[0, 0], [1, 0], [2, 0], [3, 0], [4, 0], [5, 0]]
    spacing = [1, 2, 3, 4, 5, 6]

    reduce = AlongChannelSpacingAtLeast(xy_of_node=xy_of_node, spacing=spacing)

    assert_array_equal(reduce([0, 1, 2, 3, 4, 5]), [0, 1, 3, 5])
    assert_array_equal(reduce([0, 1, 2, 3, 4, 5]), reduce(reduce([0, 1, 2, 3, 4, 5])))

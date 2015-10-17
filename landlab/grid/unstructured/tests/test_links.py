from numpy.testing import assert_array_equal
from nose.tools import raises
try:
    from nose.tools import assert_is
except ImportError:
    from landlab.testing.tools import assert_is

from landlab.grid.unstructured.links import (link_ids_at_node,
                                             in_link_ids_at_node,
                                             out_link_ids_at_node)


@raises(ValueError)
def test_array_length_mismatch():
    link_ids_at_node(([0, 1, 2, 3, 4, 5], [3, 4, 5, 6, 7]))


def test_link_arg_as_scalars():
    links, count = link_ids_at_node((0, 1), number_of_nodes=2)
    assert_array_equal(links, [0, 0])
    assert_array_equal(count, [1, 1])

    links, count = in_link_ids_at_node((0, 1), number_of_nodes=2)
    assert_array_equal(links, [0])
    assert_array_equal(count, [0, 1])

    links, count = out_link_ids_at_node((0, 1), number_of_nodes=2)
    assert_array_equal(links, [0])
    assert_array_equal(count, [1, 0])


def test_link_arg_as_pairs():
    links, count = link_ids_at_node(
        [(0, 1), (0, 2), (0, 3)], number_of_nodes=4)
    assert_array_equal(links, [0, 1, 2, 0, 1, 2])
    assert_array_equal(count, [3, 1, 1, 1])

    links, count = in_link_ids_at_node(
        [(0, 1), (0, 2), (0, 3)], number_of_nodes=4)
    assert_array_equal(links, [0, 1, 2])
    assert_array_equal(count, [0, 1, 1, 1])

    links, count = out_link_ids_at_node(
        [(0, 1), (0, 2), (0, 3)], number_of_nodes=4)
    assert_array_equal(links, [0, 1, 2])
    assert_array_equal(count, [3, 0, 0, 0])


def test_link_arg_as_tuple_sequence():
    links, count = link_ids_at_node([(0, 1), (1, 0)], number_of_nodes=2)
    assert_array_equal(links, [1, 0, 0, 1])
    assert_array_equal(count, [2, 2])

    links, count = in_link_ids_at_node([(0, 1), (1, 0)], number_of_nodes=2)
    assert_array_equal(links, [1, 0])
    assert_array_equal(count, [1, 1])

    links, count = out_link_ids_at_node([(0, 1), (1, 0)], number_of_nodes=2)
    assert_array_equal(links, [0, 1])
    assert_array_equal(count, [1, 1])


def test_links_3x3_quad():
    link_ends = [(0, 3), (1, 4), (2, 5), (3, 6), (4, 7), (5, 8),
                 (0, 1), (1, 2), (3, 4), (4, 5), (6, 7), (7, 8)]

    (links, count) = link_ids_at_node(list(zip(*link_ends)), number_of_nodes=9)

    assert_array_equal(links, [0, 6, 6, 1, 7, 7, 2,
                               0, 3, 8, 1, 8, 4, 9, 2, 9, 5,
                               3, 10, 4, 10, 11, 5, 11])
    assert_array_equal(count, [2, 3, 2, 3, 4, 3, 2, 3, 2])


def test_out_links_2x2_quad():
    link_ends = [(0, 2), (1, 3), (0, 1), (2, 3)]
    (links, count) = out_link_ids_at_node(zip(*link_ends), number_of_nodes=4)
    assert_array_equal(links, [0, 2, 1, 3])
    assert_array_equal(count, [2, 1, 1, 0])


def test_out_links_3x3_quad():
    link_ends = [(0, 3), (1, 4), (2, 5), (3, 6), (4, 7), (5, 8),
                 (0, 1), (1, 2), (3, 4), (4, 5), (6, 7), (7, 8)]

    (links, count) = out_link_ids_at_node(link_ends, number_of_nodes=9)
    assert_array_equal(links, [0, 6, 1, 7, 2, 3, 8, 4, 9, 5, 10, 11])
    assert_array_equal(count, [2, 2, 1, 2, 2, 1, 1, 1, 0])

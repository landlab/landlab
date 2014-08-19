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
    links, offset = link_ids_at_node((0, 1), number_of_nodes=2)
    assert_array_equal(links, [0, 0])
    assert_array_equal(offset, [0, 1, 2])

    links, offset = in_link_ids_at_node((0, 1), number_of_nodes=2)
    assert_array_equal(links, [0])
    assert_array_equal(offset, [0, 0, 1])

    links, offset = out_link_ids_at_node((0, 1), number_of_nodes=2)
    assert_array_equal(links, [0])
    assert_array_equal(offset, [0, 1, 1])


def test_link_arg_as_pairs():
    links, offset = link_ids_at_node([(0, 1), (0, 2), (0, 3)], number_of_nodes=4)
    assert_array_equal(links, [0, 1, 2, 0, 1, 2])
    assert_array_equal(offset, [0, 3, 4, 5, 6])

    links, offset = in_link_ids_at_node([(0, 1), (0, 2), (0, 3)], number_of_nodes=4)
    assert_array_equal(links, [0, 1, 2])
    assert_array_equal(offset, [0, 0, 1, 2, 3])

    links, offset = out_link_ids_at_node([(0, 1), (0, 2), (0, 3)], number_of_nodes=4)
    assert_array_equal(links, [0, 1, 2])
    assert_array_equal(offset, [0, 3, 3, 3, 3])


def test_link_arg_as_tuple_sequence():
    links, offset = link_ids_at_node([(0, 1), (1, 0)], number_of_nodes=2)
    assert_array_equal(links, [1, 0, 0, 1])
    assert_array_equal(offset, [0, 2, 4])

    links, offset = in_link_ids_at_node([(0, 1), (1, 0)], number_of_nodes=2)
    assert_array_equal(links, [1, 0])
    assert_array_equal(offset, [0, 1, 2])

    links, offset = out_link_ids_at_node([(0, 1), (1, 0)], number_of_nodes=2)
    assert_array_equal(links, [0, 1])
    assert_array_equal(offset, [0, 1, 2])


def test_links_3x3_quad():
    link_ends = [(0, 3), (1, 4), (2, 5), (3, 6), (4, 7), (5, 8),
                 (0, 1), (1, 2), (3, 4), (4, 5), (6, 7), (7, 8)]

    (links, offset) = link_ids_at_node(zip(*link_ends), number_of_nodes=9)
    assert_array_equal(links, [0, 6, 6, 1, 7, 7, 2,
                               0, 3, 8, 1, 8, 4, 9, 2, 9, 5,
                               3, 10, 4, 10, 11, 5, 11])
    assert_array_equal(offset, [0, 2, 5, 7, 10, 14, 17, 19, 22, 24])


def test_out_links_2x2_quad():
    link_ends = [(0, 2), (1, 3), (0, 1), (2, 3)]
    (links, offset) = out_link_ids_at_node(zip(*link_ends), number_of_nodes=4)
    assert_array_equal(links, [0, 2, 1, 3])
    assert_array_equal(offset, [0, 2, 3, 4, 4])


def test_out_links_3x3_quad():
    link_ends = [(0, 3), (1, 4), (2, 5), (3, 6), (4, 7), (5, 8),
                 (0, 1), (1, 2), (3, 4), (4, 5), (6, 7), (7, 8)]

    (links, offset) = out_link_ids_at_node(link_ends, number_of_nodes=9)
    assert_array_equal(links, [0, 6, 1, 7, 2, 3, 8, 4, 9, 5, 10, 11])
    assert_array_equal(offset, [0, 2, 4, 5, 7, 9, 10, 11, 12, 12])

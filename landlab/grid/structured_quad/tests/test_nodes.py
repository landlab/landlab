from numpy.testing import assert_array_equal

from landlab.grid.base import CLOSED_BOUNDARY, CORE_NODE, FIXED_VALUE_BOUNDARY
from landlab.grid.structured_quad import nodes
from landlab.testing.tools import assert_array_is_int


def test_perimeter_nodes():
    node_ids = nodes.perimeter((4, 5))
    assert_array_equal(node_ids, [0, 1, 2, 3, 4, 5, 9, 10, 14, 15, 16, 17, 18, 19])
    assert_array_is_int(node_ids)


def test_perimeter_status_default():
    node_status = nodes.status_with_perimeter_as_boundary((4, 5))
    F, C = CLOSED_BOUNDARY, CORE_NODE
    assert_array_equal(
        node_status,
        [[F, F, F, F, F], [F, C, C, C, F], [F, C, C, C, F], [F, F, F, F, F]],
    )
    assert_array_is_int(node_status)


def test_perimeter_status_status_as_scalar():
    node_status = nodes.status_with_perimeter_as_boundary(
        (4, 5), node_status=CLOSED_BOUNDARY
    )
    B, C = CLOSED_BOUNDARY, CORE_NODE
    assert_array_equal(
        node_status,
        [[B, B, B, B, B], [B, C, C, C, B], [B, C, C, C, B], [B, B, B, B, B]],
    )
    assert_array_is_int(node_status)


def test_perimeter_status_status_as_array():
    F, B, C = FIXED_VALUE_BOUNDARY, CLOSED_BOUNDARY, CORE_NODE
    status = [F, F, F, F, F, F, B, F, B, B, B, B, B, B]
    node_status = nodes.status_with_perimeter_as_boundary((4, 5), node_status=status)
    assert_array_equal(
        node_status,
        [[F, F, F, F, F], [F, C, C, C, B], [F, C, C, C, B], [B, B, B, B, B]],
    )
    assert_array_is_int(node_status)

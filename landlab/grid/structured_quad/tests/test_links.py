import numpy as np
import pytest
from numpy.testing import assert_array_equal

from landlab.grid.base import CLOSED_BOUNDARY, CORE_NODE
from landlab.grid.structured_quad.links import active_link_ids
from landlab.grid.structured_quad.nodes import status_with_perimeter_as_boundary
from landlab.testing.tools import assert_array_is_int


def test_active_links_ids():
    status = np.empty((4, 5), dtype=int)
    status.fill(CLOSED_BOUNDARY)
    status[1, 2] = status[1, 3] = status[2, 2] = status[2, 3] = CORE_NODE

    link_ids = active_link_ids((4, 5), status)
    assert_array_equal(link_ids, [11, 15, 16, 20])
    # assert_array_equal(link_ids, [7, 8, 21, 25])
    assert_array_is_int(link_ids)


def test_active_links_with_edge_boundaries():
    status = status_with_perimeter_as_boundary((3, 4))
    link_ids = active_link_ids((3, 4), status)
    assert_array_equal(link_ids, [8])
    assert_array_is_int(link_ids)


def test_active_link_ids_with_shape_mismatch():
    with pytest.raises(ValueError):
        active_link_ids((3, 4), np.zeros(3))

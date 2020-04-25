import numpy as np
import pytest
from numpy.testing import assert_array_equal

from landlab.grid.nodestatus import NodeStatus
from landlab.components.overland_flow._links import active_link_ids
# from landlab.testing.tools import assert_array_is_int


def test_active_links_ids():
    status = np.empty((4, 5), dtype=int)
    status.fill(NodeStatus.CLOSED)
    status[1, 2] = status[1, 3] = status[2, 2] = status[2, 3] = NodeStatus.CORE

    link_ids = active_link_ids((4, 5), status)
    assert_array_equal(link_ids, [11, 15, 16, 20])
    # assert_array_is_int(link_ids)


def test_active_link_ids_with_shape_mismatch():
    with pytest.raises(ValueError):
        active_link_ids((3, 4), np.zeros(3))

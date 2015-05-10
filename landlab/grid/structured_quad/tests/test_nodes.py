import numpy as np

from numpy.testing import assert_array_equal
from nose.tools import raises, assert_equal

from landlab.grid.structured_quad import nodes

from landlab.grid.base import CORE_NODE, FIXED_VALUE_BOUNDARY, CLOSED_BOUNDARY


def test_perimenter_nodes():
    node_ids = nodes.perimeter((4, 5))
    assert_array_equal(node_ids,
                       [0, 1, 2, 3, 4, 5, 9, 10, 14, 15, 16, 17, 18, 19])
    assert_equal(node_ids.dtype, np.int)

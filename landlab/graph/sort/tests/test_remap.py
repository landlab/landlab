from nose.tools import (assert_true, assert_false, assert_equal,
                        assert_almost_equal, assert_is, assert_is_not)
from numpy.testing import assert_array_equal, assert_array_almost_equal
import numpy as np

from landlab.graph.sort.sort import remap


def test_remap():
    src = np.array([1, 2, 3, 4])
    mapping = np.array([-1, 10, 20, 30, 40])

    rtn = remap(src, mapping)
    assert_array_equal(rtn, [10, 20, 30, 40])
    assert_is_not(rtn, src)


def test_remap_inplace():
    src = np.array([1, 2, 3, 4])
    mapping = np.array([-1, 10, 20, 30, 40])

    rtn = remap(src, mapping, inplace=True)

    assert_array_equal(rtn, [10, 20, 30, 40])
    assert_is(rtn, src)


def test_remap_out():
    src = np.array([1, 2, 3, 4])
    dst = np.empty_like(src)
    mapping = np.array([-1, 10, 20, 30, 40])

    rtn = remap(src, mapping, out=dst)

    assert_array_equal(rtn, [10, 20, 30, 40])
    assert_is(rtn, dst)

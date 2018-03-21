"""Test the flow accumulator component.

@author: krb
"""
# Created on Thurs Nov 12, 2015
import os

import numpy as np
from numpy.testing import assert_array_equal, assert_array_almost_equal

from nose.tools import (with_setup, assert_true, assert_false, assert_raises,
                        assert_almost_equal, assert_equal)
try:
    from nose.tools import (assert_is, assert_set_equal, assert_dict_equal)
except ImportError:
    from landlab.testing.tools import (assert_is, assert_set_equal,
                                       assert_dict_equal)

import landlab
from landlab import RasterModelGrid, HexModelGrid, FieldError # NetworkModelGrid,

from landlab import CLOSED_BOUNDARY
from landlab import BAD_INDEX_VALUE as XX

# def test_NetworkModelGrid():
#     """ """
#     pass


def test_RasterModelGrid():
    """ """
    pass


def test_HexModelGrid():
    """ """
    pass


def test_different_length_data():
    """ """
    pass


def test_bad_location_index():
    """ """
    pass


def test_bad_grid_element_name():
    """ """
    pass


def test_grid_element_size():
    """ """
    pass

def test_different_functions():
    """ """
    pass
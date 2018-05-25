#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Apr 30 21:28:33 2018

@author: barnhark
"""

#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Apr 30 09:17:36 2018

@author: barnhark
"""

import numpy as np
#from numpy.testing import assert_array_equal, assert_array_almost_equal
from nose.tools import assert_raises#, assert_almost_equal, assert_equal

from landlab import RasterModelGrid
from landlab.components import LithoLayers

def test_z0s_ids_different_shape():
    """Test that providing z0s and ids of different shapes raises an error."""
    mg = RasterModelGrid(3, 3)
    z = mg.add_zeros('node', 'topographic__elevation')
    z0s = [-4, -3, -2, -1, 0, 1, 2, 3, 4]
    ids = [1, 2, 1, 2, 1, 2, 1, 2, 1, 2]
    attrs = {'K_sp': {1: 0.001, 2: 0.0001}}
    assert_raises(ValueError, LithoLayers, mg, z0s, ids, attrs)


def test_z0s_bad_order():
    """Test that providing z0s in a bad order raises an error."""
    mg = RasterModelGrid(3, 3)
    z = mg.add_zeros('node', 'topographic__elevation')
    z0s = [-4, -3, -2, -1, 0, 1, 2, 6, 4]
    ids = [1, 2, 1, 2, 1, 2, 1, 2, 1]
    attrs = {'K_sp': {1: 0.001, 2: 0.0001}}
    assert_raises(ValueError, LithoLayers, mg, z0s, ids, attrs)


def test_bad_function():
    """Test that providing a function of three variables."""
    mg = RasterModelGrid(3, 3)
    z = mg.add_zeros('node', 'topographic__elevation')
    z0s = [-4, -3, -2, -1, 0, 1, 2, 4, 6]
    ids = [1, 2, 1, 2, 1, 2, 1, 2, 1]
    attrs = {'K_sp': {1: 0.001, 2: 0.0001}}
    func = lambda x, y, z: 0*x + 0*y + z
    assert_raises(ValueError, LithoLayers, mg, z0s, ids, attrs, function=func)

def test_function_returns_wrong_number_of_values():
    """Test that providing a function that returns one value raises error."""
    mg = RasterModelGrid(3, 3)
    z = mg.add_zeros('node', 'topographic__elevation')
    z0s = [-4, -3, -2, -1, 0, 1, 2, 4, 6]
    ids = [1, 2, 1, 2, 1, 2, 1, 2, 1]
    attrs = {'K_sp': {1: 0.001, 2: 0.0001}}
    func = lambda x, y, z: np.mean(x+y)
    assert_raises(ValueError, LithoLayers, mg, z0s, ids, attrs, function=func)

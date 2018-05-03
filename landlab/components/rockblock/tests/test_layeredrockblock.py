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

#import numpy as np
#from numpy.testing import assert_array_equal, assert_array_almost_equal
from nose.tools import assert_raises#, assert_almost_equal, assert_equal

from landlab import RasterModelGrid
from landlab.components import LayeredRockBlock

def test_z0s_ids_different_shape():
    """Test that providing z0s and ids of different shapes raises an error."""
    mg = RasterModelGrid(3, 3)
    z = mg.add_zeros('node', 'topographic__elevation')
    z0s = [-4, -3, -2, -1, 0, 1, 2, 3, 4]
    ids = [1, 2, 1, 2, 1, 2, 1, 2, 1, 2]
    attrs = {'K_sp': {1: 0.001, 2: 0.0001}}
    assert_raises(ValueError, LayeredRockBlock, mg, z0s, ids, attrs)


def test_z0s_bad_order():
    """Test that providing z0s in a bad order raises an error."""
    mg = RasterModelGrid(3, 3)
    z = mg.add_zeros('node', 'topographic__elevation')
    z0s = [-4, -3, -2, -1, 0, 1, 2, 6, 4]
    ids = [1, 2, 1, 2, 1, 2, 1, 2, 1]
    attrs = {'K_sp': {1: 0.001, 2: 0.0001}}
    assert_raises(ValueError, LayeredRockBlock, mg, z0s, ids, attrs)

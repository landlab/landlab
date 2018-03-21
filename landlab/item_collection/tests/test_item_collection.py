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
from landlab.item_collection import ItemCollection

from landlab import CLOSED_BOUNDARY
from landlab import BAD_INDEX_VALUE as XX

_LOCATIONS = ['node', 'patch', 'link', 'corner', 'face', 'cell']

def test_RasterAndHexModelGrids():
    """Test compatability with RasterModelGrid and HexModelGrid"""

    element_id = [0, 0, 1, 1, 2, 3, 5]
    volume = [1, 2, 3, 4, 5, 6, 7]
    age = [10, 11, 12, 13, 14, 15, 16] 
    
    rgrid = RasterModelGrid(6,6)
    hgrid = HexModelGrid(4,5)
    # ngrid = NetworkModelGrid(...)
    
    # for each grid type
    for grid in [rgrid, hgrid]: # [rgrid, hgrid, ngrid]
        # for each location, test that you can make an item collection there:
        
#        if isinstance(grid, NetworkModelGrid):
#            locations = ['node', 'link']
#        else:
#            locations = _LOCATIONS
            
        locations = _LOCATIONS
        for loc in locations:
            # recreate data
            data = {'age': age,
                    'volume': volume}
        
            # test that an item collection can be created
            ic = ItemCollection(grid, 
                                data = data, 
                                grid_element = loc, 
                                element_id = element_id)
        
            # test that a bad grid element name raises the right error
            data = {'age': age,
                    'volume': volume}
            assert_raises(ValueError, 
                          ItemCollection, 
                          grid, 
                          data = data, 
                          grid_element = 'foo', 
                          element_id = element_id)
            
            # test that an too big element raises a warning
            data = {'age': age,
                    'volume': volume}
            bad_element_id = [0, 0, 1, 1, 2, 3, 500]
            assert_raises(ValueError, 
                          ItemCollection, 
                          grid, 
                          data = data, 
                          grid_element = loc, 
                          element_id = bad_element_id)
            
            # test that a too small element raises a warming
            data = {'age': age,
                    'volume': volume}
            bad_element_id = [0, 0, 1, 1, 2, 3, -1]
            assert_raises(ValueError, 
                          ItemCollection, 
                          grid, 
                          data = data, 
                          grid_element = loc, 
                          element_id = bad_element_id)      
            
            # test wrong length data
            data = {'age': age[:-1],
                    'volume': volume}
            assert_raises(ValueError, 
                          ItemCollection, 
                          grid, 
                          data = data, 
                          grid_element = loc, 
                          element_id = bad_element_id) 
            

def test_grid_element_size():
    """ """
    pass

def test_different_functions():
    """ """
    pass
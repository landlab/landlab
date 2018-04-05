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


element_id = [0, 0, 1, 1, 2, 3, 5]
volume = [1, 2, 3, 4, 5, 6, 7]
age = [10, 11, 12, 13, 14, 15, 16]

rgrid = RasterModelGrid(6,6)
hgrid = HexModelGrid(4,5)
# ngrid = NetworkModelGrid(...)

def test_instantiation():
    """Test instantiation"""
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

            # test that grid_element can be a list
            size = len(volume)
            grid_element = np.empty(size, dtype=object)
            grid_element.fill(loc)
            data = {'age': age,
                    'volume': volume}

            ic = ItemCollection(grid,
                                data = data,
                                grid_element = grid_element,
                                element_id = element_id)

        # test that grid element can be mutiple types
        size = len(volume)
        grid_element = np.empty(size, dtype=object)
        grid_element.fill('node')
        grid_element[0] = 'link'
        data = {'age': age,
                'volume': volume}

        ic = ItemCollection(grid,
                            data = data,
                            grid_element = grid_element,
                            element_id = element_id)

def test_bad_grid_element():
    """Test bad grid element."""
    # for each grid type
    for grid in [rgrid, hgrid]: # [rgrid, hgrid, ngrid]
        # for each location, test that you can make an item collection there:

#        if isinstance(grid, NetworkModelGrid):
#            locations = ['node', 'link']
#        else:
#            locations = _LOCATIONS

        locations = _LOCATIONS
        for loc in locations:
            data = {'age': age,
                    'volume': volume}

            assert_raises(ValueError,
                          ItemCollection,
                          grid,
                          data = data,
                          grid_element = 'foo',
                          element_id = element_id)

            size = len(volume)
            grid_element = np.empty(size, dtype=object)
            grid_element.fill('node')
            grid_element[0] = 'foo'
            data = {'age': age,
                    'volume': volume}

            assert_raises(ValueError,
                          ItemCollection,
                          grid,
                          data = data,
                          grid_element = grid_element,
                          element_id = element_id)

def test_non_int_element_id():
    """Test a grid element that is non-integer."""
    bad_element_id = [0.0, 0.0, 1.1, 1.2, 2.3, 3.0, 500.0]
    # for each grid type
    for grid in [rgrid, hgrid]: # [rgrid, hgrid, ngrid]
        # for each location, test that you can make an item collection there:

#        if isinstance(grid, NetworkModelGrid):
#            locations = ['node', 'link']
#        else:
#            locations = _LOCATIONS

        locations = _LOCATIONS
        for loc in locations:
            # test that an too big element raises a warning
            data = {'age': age,
                    'volume': volume}

            assert_raises(ValueError,
                          ItemCollection,
                          grid,
                          data = data,
                          grid_element = loc,
                          element_id = bad_element_id)

def test_big_element_id():
    """Test a grid element that is too big."""
    bad_element_id = [0, 0, 1, 1, 2, 3, 500]
    # for each grid type
    for grid in [rgrid, hgrid]: # [rgrid, hgrid, ngrid]
        # for each location, test that you can make an item collection there:

#        if isinstance(grid, NetworkModelGrid):
#            locations = ['node', 'link']
#        else:
#            locations = _LOCATIONS

        locations = _LOCATIONS
        for loc in locations:
            # test that an too big element raises a warning
            data = {'age': age,
                    'volume': volume}

            assert_raises(ValueError,
                          ItemCollection,
                          grid,
                          data = data,
                          grid_element = loc,
                          element_id = bad_element_id)


def test_small_element_id():
    """Test a grid element that is too small."""
    bad_element_id = [0, 0, 1, 1, 2, 3, -1]
    # for each grid type
    for grid in [rgrid, hgrid]: # [rgrid, hgrid, ngrid]
        # for each location, test that you can make an item collection there:

#        if isinstance(grid, NetworkModelGrid):
#            locations = ['node', 'link']
#        else:
#            locations = _LOCATIONS

        locations = _LOCATIONS
        for loc in locations:
            # test that a too small element raises a warming
            data = {'age': age,
                    'volume': volume}
            assert_raises(ValueError,
                          ItemCollection,
                          grid,
                          data = data,
                          grid_element = loc,
                          element_id = bad_element_id)

def test_wrong_length_data():
    """Test a grid element that is too small."""
    # for each grid type
    for grid in [rgrid, hgrid]: # [rgrid, hgrid, ngrid]
        # for each location, test that you can make an item collection there:

#        if isinstance(grid, NetworkModelGrid):
#            locations = ['node', 'link']
#        else:
#            locations = _LOCATIONS

        locations = _LOCATIONS
        for loc in locations:
            # test that a too small element raises a warming
            data = {'age': age[:-1],
                    'volume': volume}
            # test wrong length data
            assert_raises(ValueError,
                          ItemCollection,
                          grid,
                          data = data,
                          grid_element = loc,
                          element_id = element_id)


def test_element_id_size():
    """Test passing the wrong size element_id"""
    # for each grid type
    for grid in [rgrid, hgrid]: # [rgrid, hgrid, ngrid]
        # for each location, test that you can make an item collection there:

#        if isinstance(grid, NetworkModelGrid):
#            locations = ['node', 'link']
#        else:
#            locations = _LOCATIONS

        locations = _LOCATIONS
        for loc in locations:
            # test that a too small element raises a warming
            data = {'age': age[:-1],
                    'volume': volume}

            size = len(volume)
            grid_element = np.empty(size, dtype=object)
            grid_element.fill(loc)

            # test wrong length data
            assert_raises(ValueError,
                          ItemCollection,
                          grid,
                          data = data,
                          grid_element = grid_element,
                          element_id = element_id[:-1])

def test_grid_element_size():
    """Test passing the wrong size grid element"""
    # for each grid type
    for grid in [rgrid, hgrid]: # [rgrid, hgrid, ngrid]
        # for each location, test that you can make an item collection there:

#        if isinstance(grid, NetworkModelGrid):
#            locations = ['node', 'link']
#        else:
#            locations = _LOCATIONS

        locations = _LOCATIONS
        for loc in locations:
            # test that a too small element raises a warming
            data = {'age': age[:-1],
                    'volume': volume}

            size = len(volume)
            grid_element = np.empty(size, dtype=object)
            grid_element.fill(loc)

            # test wrong length data
            assert_raises(ValueError,
                          ItemCollection,
                          grid,
                          data = data,
                          grid_element = grid_element[:-1],
                          element_id = element_id)


def test_adding_bad_size_variable():
    """Test passing the wrong size variable."""
    new_var = [1, 2, 3, 4, 5, 6]

    # for each grid type
    for grid in [rgrid, hgrid]: # [rgrid, hgrid, ngrid]
        # for each location, test that you can make an item collection there:

#        if isinstance(grid, NetworkModelGrid):
#            locations = ['node', 'link']
#        else:
#            locations = _LOCATIONS

        locations = _LOCATIONS
        for loc in locations:
            # test that a too small element raises a warming
            data = {'age': age,
                    'volume': volume}


            # test wrong length data
            ic = ItemCollection(grid,
                                data = data,
                                grid_element = loc,
                                element_id = element_id)

            assert_raises(ValueError,
                          ic.add_variable,
                          'new_var',
                          new_var)


def test_adding_non_string_variable_name():
    """Test passing a variable name that is not a string."""
    new_var = [10, 11, 12, 13, 14, 15, 16]

    # for each grid type
    for grid in [rgrid, hgrid]: # [rgrid, hgrid, ngrid]
        # for each location, test that you can make an item collection there:

#        if isinstance(grid, NetworkModelGrid):
#            locations = ['node', 'link']
#        else:
#            locations = _LOCATIONS

        locations = _LOCATIONS
        for loc in locations:
            # test that a too small element raises a warming
            data = {'age': age,
                    'volume': volume}

            size = len(volume)
            grid_element = np.empty(size, dtype=object)
            grid_element.fill(loc)

            # test wrong length data
            ic = ItemCollection(grid,
                                data = data,
                                grid_element = loc,
                                element_id = element_id)

            assert_raises(ValueError,
                          ic.add_variable,
                          1.1,
                          new_var)


def test_adding_old_variable_name():
    """Test adding an old variable name."""
    new_var = [10, 11, 12, 13, 14, 15, 16]

    # for each grid type
    for grid in [rgrid, hgrid]: # [rgrid, hgrid, ngrid]
        # for each location, test that you can make an item collection there:

#        if isinstance(grid, NetworkModelGrid):
#            locations = ['node', 'link']
#        else:
#            locations = _LOCATIONS

        locations = _LOCATIONS
        for loc in locations:
            # test that a too small element raises a warming
            data = {'age': age,
                    'volume': volume}

            # test wrong length data
            ic = ItemCollection(grid,
                                data = data,
                                grid_element = loc,
                                element_id = element_id)

            assert_raises(ValueError,
                          ic.add_variable,
                          'age',
                          new_var)


def test_adding_items_with_extra_variables():
    """Test adding items with extra variables."""
    new_age = [0, 1]
    new_volume = [0, 3]
    foo = [4, 5]
    # for each grid type
    for grid in [rgrid, hgrid]: # [rgrid, hgrid, ngrid]
        # for each location, test that you can make an item collection there:

#        if isinstance(grid, NetworkModelGrid):
#            locations = ['node', 'link']
#        else:
#            locations = _LOCATIONS

        locations = _LOCATIONS
        for loc in locations:
            # test that a too small element raises a warming
            data = {'age': age,
                    'volume': volume}

            # test wrong length data
            ic = ItemCollection(grid,
                                data = data,
                                grid_element = loc,
                                element_id = element_id)
            new_data = {'age': new_age,
                        'volume': new_volume,
                        'foo': foo}

            assert_raises(ValueError,
                          ic.add_item,
                          data = new_data,
                          grid_element = loc,
                          element_id = [6, 7])


def test_adding_items_with_not_enough_variables():
    """Test adding items with not_enough variables."""
    new_age = [0, 1]
    # for each grid type
    for grid in [rgrid, hgrid]: # [rgrid, hgrid, ngrid]
        # for each location, test that you can make an item collection there:

#        if isinstance(grid, NetworkModelGrid):
#            locations = ['node', 'link']
#        else:
#            locations = _LOCATIONS

        locations = _LOCATIONS
        for loc in locations:
            # test that a too small element raises a warming
            data = {'age': age,
                    'volume': volume}

            # test wrong length data
            ic = ItemCollection(grid,
                                data = data,
                                grid_element = loc,
                                element_id = element_id)
            new_data = {'age': new_age}

            assert_raises(ValueError,
                          ic.add_item,
                          data = new_data,
                          grid_element = loc,
                          element_id = [6, 7])

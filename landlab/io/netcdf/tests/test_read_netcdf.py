#! /usr/bin/env python
"""Unit tests for landlab.io.netcdf module."""

import os
from nose.tools import assert_equal
from nose.plugins.skip import SkipTest

from landlab.io.netcdf import read_netcdf, WITH_NETCDF4


_TEST_DATA_DIR = os.path.join(os.path.dirname(__file__), 'data')

grid_mapping_keys = ['grid_mapping_name', 'longitude_of_central_meridian', 
                     'false_easting', 'false_northing', 
                     'latitude_of_projection_origin', 
                     'scale_factor_at_central_meridian', 'long_name', 
                     'longitude_of_prime_meridian', 'semi_major_axis', 
                     'inverse_flattening', 'spatial_ref', 'GeoTransform']

def test_read_netcdf3_64bit():
    """Test read_netcdf for with 64-bit netcdf3 format."""
    grid = read_netcdf(os.path.join(_TEST_DATA_DIR, 'test-netcdf3-64bit.nc'))
    assert_equal(grid.shape, (4, 3))


def test_read_netcdf4():
    """Test read_netcdf with netcdf4 format."""
    if not WITH_NETCDF4:
        raise SkipTest('netCDF4 package not installed')

    grid = read_netcdf(os.path.join(_TEST_DATA_DIR, 'test-netcdf4.nc'))
    assert_equal(grid.shape, (4, 3))


def test_netcdf_with_grid_mapping_3():
    """Test read netcdf with grid mapping.""" 
    grid = read_netcdf(os.path.join(_TEST_DATA_DIR, 'grid_mapping_ex.nc'))
    assert_equal(hasattr(grid, 'grid_mapping'), True)
    mapping = grid.grid_mapping
    assert_equal(type(mapping), dict)
    for gmk in grid_mapping_keys:
        assert_equal(gmk in mapping, True)
    
    
def test_netcdf_with_grid_mapping_4():
    """Test read netcdf with grid mapping.""" 
    grid = read_netcdf(os.path.join(_TEST_DATA_DIR, 'grid_mapping_ex_NETCDF4.nc'))
    assert_equal(hasattr(grid, 'grid_mapping'), True)
    mapping = grid.grid_mapping
    assert_equal(type(mapping), dict)
    for gmk in grid_mapping_keys:
        assert_equal(gmk in mapping, True)
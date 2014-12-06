#! /usr/bin/env python
"""
Unit tests for landlab.io.netcdf module.
"""

import os
import numpy as np
from StringIO import StringIO
from unittest import skipIf
from nose.tools import assert_equal, assert_true, assert_raises
try:
    from nose import assert_list_equal
except ImportError:
    from landlab.testing.tools import assert_list_equal

from landlab import RasterModelGrid
from landlab.io.netcdf import write_netcdf, NotRasterGridError, WITH_NETCDF4
from landlab.io.netcdf.read import _get_raster_spacing


_TEST_DATA_DIR = os.path.join(os.path.dirname(__file__), 'data')


@skipIf(not WITH_NETCDF4, 'netCDF4 package not installed')
def test_netcdf_write():
    field = RasterModelGrid(4, 3)
    #field.new_field_location('node', 12.)
    values = np.arange(12.)
    field.add_field('node', 'topographic_elevation', values)

    write_netcdf('test.nc', field)

    import netCDF4 as nc

    root = nc.Dataset('test.nc', 'r', format='NETCDF4')
    assert_equal(set(root.dimensions), set(['ni', 'nj', 'nt']))
    assert_equal(len(root.dimensions['ni']), 3)
    assert_equal(len(root.dimensions['nj']), 4)
    assert_true(len(root.dimensions['nt']), 1)
    assert_true(root.dimensions['nt'].isunlimited())

    assert_equal(set(root.variables),
                 set(['x', 'y', 'topographic_elevation']))

    assert_list_equal(list(root.variables['x'][:].flat),
                      [0., 1., 2.,
                       0., 1., 2.,
                       0., 1., 2.,
                       0., 1., 2., ])
    assert_list_equal(list(root.variables['y'][:].flat),
                      [0., 0., 0.,
                       1., 1., 1.,
                       2., 2., 2.,
                       3., 3., 3., ])
    assert_list_equal(
        list(root.variables['topographic_elevation'][:].flat),
        range(12))
    root.close()


def test_2d_unit_spacing():
    (x, y) = np.meshgrid(np.arange(5.), np.arange(4.))

    spacing = _get_raster_spacing((y, x))
    assert_equal(spacing, 1.)

def test_2d_non_unit_spacing():
    (x, y) = np.meshgrid(np.arange(5.) * 2, np.arange(4.) * 2)

    spacing = _get_raster_spacing((y, x))
    assert_equal(spacing, 2.)

def test_2d_uneven_spacing_axis_0():
    (x, y) = np.meshgrid(np.logspace(0., 2., num=5), np.arange(4.))

    assert_raises(NotRasterGridError, _get_raster_spacing, (y, x))

def test_2d_uneven_spacing_axis_1():
    (x, y) = np.meshgrid(np.arange(4.), np.logspace(0., 2., num=5))

    assert_raises(NotRasterGridError, _get_raster_spacing, (y, x))

def test_2d_switched_coords():
    (x, y) = np.meshgrid(np.arange(5.), np.arange(4.))

    spacing = _get_raster_spacing((x, y))
    assert_equal(spacing, 0.)


def test_1d__unit_spacing():
    spacing = _get_raster_spacing((np.arange(5.), ))
    assert_equal(spacing, 1.)

def test_1d_non_unit_spacing():
    spacing = _get_raster_spacing((np.arange(5.) * 2, ))
    assert_equal(spacing, 2.)

def test_1d_uneven_spacing():
    (x, y) = np.meshgrid(np.logspace(0., 2., num=5), np.arange(4.))

    assert_raises(NotRasterGridError, _get_raster_spacing, 
                  (np.logspace(0., 2., num=5), ))

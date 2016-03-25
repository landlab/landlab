#! /usr/bin/env python
import os

import numpy as np
from numpy.testing import assert_array_almost_equal
from nose.tools import assert_true, assert_equal, assert_raises
try:
    from nose.tools import assert_list_equal
except ImportError:
    from landlab.testing.tools import assert_list_equal

from landlab.testing.tools import cdtemp
from landlab.io import read_esri_ascii
from landlab.io.netcdf import read_netcdf
from landlab import RasterModelGrid


def test_save_esri_ascii():
    grid = RasterModelGrid(4, 5, 2.)
    grid.add_field('node', 'air__temperature', np.arange(20.))

    with cdtemp() as _:
        grid.save('test.asc', format='esri-ascii')
        assert_true(os.path.isfile('test.asc'))


def test_add_extension():
    grid = RasterModelGrid(4, 5, 2.)
    grid.add_field('node', 'air__temperature', np.arange(20.))

    with cdtemp() as _:
        grid.save('test', format='esri-ascii')
        assert_true(os.path.isfile('test.asc'))

    with cdtemp() as _:
        grid.save('test', format='netcdf')
        assert_true(os.path.isfile('test.nc'))


def test_replace_extension():
    grid = RasterModelGrid(4, 5, 2.)
    grid.add_field('node', 'air__temperature', np.arange(20.))

    with cdtemp() as _:
        grid.save('test.nc', format='esri-ascii')
        assert_true(os.path.isfile('test.asc'))

    with cdtemp() as _:
        grid.save('test.asc', format='netcdf')
        assert_true(os.path.isfile('test.nc'))


def test_guess_format():
    grid = RasterModelGrid(4, 5, 2.)
    grid.add_field('node', 'air__temperature', np.arange(20.))

    with cdtemp() as _:
        grid.save('test.asc')
        assert_true(os.path.isfile('test.asc'))
        read_esri_ascii('test.asc')

    with cdtemp() as _:
        grid.save('test.nc')
        assert_true(os.path.isfile('test.nc'))
        read_netcdf('test.nc')


def test_names_keyword():
    grid = RasterModelGrid(4, 5, 2.)
    grid.add_field('node', 'air__temperature', np.arange(20.))
    grid.add_field('node', 'land_surface__elevation', np.arange(20.))

    with cdtemp() as _:
        grid.save('test.asc', names='land_surface__elevation')
        assert_true(os.path.isfile('test.asc'))
        read_esri_ascii('test.asc')

    with cdtemp() as _:
        grid.save('test.asc', names=['land_surface__elevation',
                                     'air__temperature'])
        files = ['test_land_surface__elevation.asc',
                 'test_air__temperature.asc']
        for fname in files:
            assert_true(os.path.isfile(fname))
            read_esri_ascii(fname)

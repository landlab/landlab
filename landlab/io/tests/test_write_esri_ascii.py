#! /usr/bin/env python
import os

import numpy as np

from nose.tools import assert_true, assert_raises
try:
    from nose.tools import assert_list_equal
except ImportError:
    from landlab.testing.tools import assert_list_equal

from landlab.testing.tools import cdtemp
from landlab.io import write_esri_ascii
from landlab import RasterModelGrid


def test_grid_with_no_fields():
    grid = RasterModelGrid(4, 5, dx=2.)
    with cdtemp() as _:
        with assert_raises(ValueError):
            write_esri_ascii('test.asc', grid)


def test_grid_with_one_field():
    grid = RasterModelGrid(4, 5, dx=2.)
    grid.add_field('node', 'air__temperature', np.arange(20.))
    with cdtemp() as _:
        files = write_esri_ascii('test.asc', grid)
        assert_list_equal(files, ['test.asc'])
        for fname in files:
            assert_true(os.path.isfile(fname))


def test_grid_with_two_fields():
    grid = RasterModelGrid(4, 5, dx=2.)
    grid.add_field('node', 'air__temperature', np.arange(20.))
    grid.add_field('node', 'land_surface__elevation', np.arange(20.))
    with cdtemp() as _:
        files = write_esri_ascii('test.asc', grid)
        files.sort()
        assert_list_equal(files, ['test_air__temperature.asc',
                                  'test_land_surface__elevation.asc'])
        for fname in files:
            assert_true(os.path.isfile(fname))


def test_names_keyword_as_str_or_list():
    grid = RasterModelGrid(4, 5, dx=2.)
    grid.add_field('node', 'air__temperature', np.arange(20.))
    grid.add_field('node', 'land_surface__elevation', np.arange(20.))

    with cdtemp() as _:
        files = write_esri_ascii('test.asc', grid, names='air__temperature')
        assert_list_equal(files, ['test.asc'])
        assert_true(os.path.isfile('test.asc'))

    with cdtemp() as _:
        files = write_esri_ascii('test.asc', grid, names=['air__temperature'])
        assert_list_equal(files, ['test.asc'])
        assert_true(os.path.isfile('test.asc'))


def test_names_keyword_multiple_names():
    grid = RasterModelGrid(4, 5, dx=2.)
    grid.add_field('node', 'air__temperature', np.arange(20.))
    grid.add_field('node', 'land_surface__elevation', np.arange(20.))

    with cdtemp() as _:
        files = write_esri_ascii('test.asc', grid,
                                 names=['air__temperature',
                                        'land_surface__elevation'])
        files.sort()
        assert_list_equal(files, ['test_air__temperature.asc',
                                  'test_land_surface__elevation.asc'])
        for fname in files:
            assert_true(os.path.isfile(fname))


def test_names_keyword_with_bad_name():
    grid = RasterModelGrid(4, 5, dx=2.)
    grid.add_field('node', 'air__temperature', np.arange(20.))

    with cdtemp() as _:
        with assert_raises(ValueError):
            write_esri_ascii('test.asc', grid, names='not_a_name')


def test_clobber_keyword():
    grid = RasterModelGrid(4, 5, dx=2.)
    grid.add_field('node', 'air__temperature', np.arange(20.))

    with cdtemp() as _:
        write_esri_ascii('test.asc', grid)
        with assert_raises(ValueError):
            write_esri_ascii('test.asc', grid)
        with assert_raises(ValueError):
            write_esri_ascii('test.asc', grid, clobber=False)
        write_esri_ascii('test.asc', grid, clobber=True)

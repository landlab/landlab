#! /usr/bin/env python
import os

import pytest
import numpy as np
from numpy.testing import assert_array_almost_equal

from landlab.testing.tools import cdtemp
from landlab.io import write_esri_ascii, read_esri_ascii
from landlab import RasterModelGrid


def test_grid_with_no_fields():
    grid = RasterModelGrid((4, 5), spacing=(2., 2.))
    with cdtemp() as _:
        with pytest.raises(ValueError):
            write_esri_ascii('test.asc', grid)


def test_grid_with_one_field():
    grid = RasterModelGrid((4, 5), spacing=(2., 2.))
    grid.add_field('node', 'air__temperature', np.arange(20.))
    with cdtemp() as _:
        files = write_esri_ascii('test.asc', grid)
        assert files == ['test.asc']
        for fname in files:
            assert os.path.isfile(fname)


def test_grid_with_two_fields():
    grid = RasterModelGrid((4, 5), spacing=(2., 2.))
    grid.add_field('node', 'air__temperature', np.arange(20.))
    grid.add_field('node', 'land_surface__elevation', np.arange(20.))
    with cdtemp() as _:
        files = write_esri_ascii('test.asc', grid)
        files.sort()
        assert files == ['test_air__temperature.asc',
                         'test_land_surface__elevation.asc']
        for fname in files:
            assert os.path.isfile(fname)


def test_names_keyword_as_str_or_list():
    grid = RasterModelGrid((4, 5), spacing=(2., 2.))
    grid.add_field('node', 'air__temperature', np.arange(20.))
    grid.add_field('node', 'land_surface__elevation', np.arange(20.))

    with cdtemp() as _:
        files = write_esri_ascii('test.asc', grid, names='air__temperature')
        assert files == ['test.asc']
        assert os.path.isfile('test.asc')

    with cdtemp() as _:
        files = write_esri_ascii('test.asc', grid, names=['air__temperature'])
        assert files == ['test.asc']
        assert os.path.isfile('test.asc')


def test_names_keyword_multiple_names():
    grid = RasterModelGrid((4, 5), spacing=(2., 2.))
    grid.add_field('node', 'air__temperature', np.arange(20.))
    grid.add_field('node', 'land_surface__elevation', np.arange(20.))

    with cdtemp() as _:
        files = write_esri_ascii('test.asc', grid,
                                 names=['air__temperature',
                                        'land_surface__elevation'])
        files.sort()
        assert files == ['test_air__temperature.asc',
                         'test_land_surface__elevation.asc']
        for fname in files:
            assert os.path.isfile(fname)


def test_names_keyword_with_bad_name():
    grid = RasterModelGrid((4, 5), spacing=(2., 2.))
    grid.add_field('node', 'air__temperature', np.arange(20.))

    with cdtemp() as _:
        with pytest.raises(ValueError):
            write_esri_ascii('test.asc', grid, names='not_a_name')


def test_clobber_keyword():
    grid = RasterModelGrid((4, 5), spacing=(2., 2.))
    grid.add_field('node', 'air__temperature', np.arange(20.))

    with cdtemp() as _:
        write_esri_ascii('test.asc', grid)
        with pytest.raises(ValueError):
            write_esri_ascii('test.asc', grid)
        with pytest.raises(ValueError):
            write_esri_ascii('test.asc', grid, clobber=False)
        write_esri_ascii('test.asc', grid, clobber=True)


def test_write_then_read():
    grid = RasterModelGrid((4, 5), spacing=(2., 2.))
    grid.add_field('node', 'air__temperature', np.arange(20.))

    with cdtemp() as _:
        write_esri_ascii('test.asc', grid)
        new_grid, field = read_esri_ascii('test.asc')

    assert grid.number_of_node_columns == new_grid.number_of_node_columns
    assert grid.number_of_node_rows == new_grid.number_of_node_rows
    assert grid.dx == new_grid.dx
    assert_array_almost_equal(grid.node_x, new_grid.node_x)
    assert_array_almost_equal(grid.node_y, new_grid.node_y)
    assert_array_almost_equal(field, grid.at_node['air__temperature'])

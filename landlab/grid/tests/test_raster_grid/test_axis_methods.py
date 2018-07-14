import pytest
import numpy as np
from numpy.testing import assert_array_equal
from nose.tools import with_setup

from landlab import RasterModelGrid


def setup_default_grid():
    globals().update({
        'rmg': RasterModelGrid(4, 5)
    })


def setup_lon_lat_grid():
    globals().update({
        'rmg': RasterModelGrid(4, 5,
                               axis_name=['longitude', 'latitude'],
                               axis_units=['degrees_east', 'degrees_north'])
    })


@with_setup(setup_default_grid)
def test_default_names():
    assert rmg.axis_name == ('y', 'x')


@with_setup(setup_lon_lat_grid)
def test_name_keyword():
    assert rmg.axis_name == ('longitude', 'latitude')


@with_setup(setup_lon_lat_grid)
def test_name_setter():
    rmg.axis_name = ('yyy', 'xxx')
    assert rmg.axis_name == ('yyy', 'xxx')


@with_setup(setup_default_grid)
def test_name_setter_too_few_names():
    with pytest.raises(ValueError):
        rmg.axis_name = ('z', )


@with_setup(setup_default_grid)
def test_name_setter_too_many_names():
    with pytest.raises(ValueError):
        rmg.axis_name = ('z', 'y', 'x')


@with_setup(setup_default_grid)
def test_default_units():
    assert rmg.axis_units == ('-', '-')


@with_setup(setup_lon_lat_grid)
def test_name_keyword():
    assert rmg.axis_units == ('degrees_east', 'degrees_north')


@with_setup(setup_lon_lat_grid)
def test_name_setter():
    rmg.axis_units = ('mm', 'cm')
    assert rmg.axis_units == ('mm', 'cm')


@with_setup(setup_default_grid)
def test_name_setter_too_few_units():
    with pytest.raises(ValueError):
        rmg.axis_units = ('m', )


@with_setup(setup_default_grid)
def test_name_setter_too_many_units():
    with pytest.raises(ValueError):
        rmg.axis_units = ('m', 'cm', 'mm')

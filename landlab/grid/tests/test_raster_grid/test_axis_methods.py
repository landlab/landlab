import numpy as np
from numpy.testing import assert_array_equal
from nose.tools import with_setup, raises
try:
    from nose.tools import assert_tuple_equal
except ImportError:
    from landlab.testing.tools import assert_tuple_equal


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
    assert_tuple_equal(rmg.axis_name, ('y', 'x'))


@with_setup(setup_lon_lat_grid)
def test_name_keyword():
    assert_tuple_equal(rmg.axis_name, ('longitude', 'latitude'))


@with_setup(setup_lon_lat_grid)
def test_name_setter():
    rmg.axis_name = ('yyy', 'xxx')
    assert_tuple_equal(rmg.axis_name, ('yyy', 'xxx'))


@raises(ValueError)
@with_setup(setup_default_grid)
def test_name_setter_too_few_names():
    rmg.axis_name = ('z', )


@raises(ValueError)
@with_setup(setup_default_grid)
def test_name_setter_too_many_names():
    rmg.axis_name = ('z', 'y', 'x')


@with_setup(setup_default_grid)
def test_default_units():
    assert_tuple_equal(rmg.axis_units, ('-', '-'))


@with_setup(setup_lon_lat_grid)
def test_name_keyword():
    assert_tuple_equal(rmg.axis_units, ('degrees_east', 'degrees_north'))


@with_setup(setup_lon_lat_grid)
def test_name_setter():
    rmg.axis_units = ('mm', 'cm')
    assert_tuple_equal(rmg.axis_units, ('mm', 'cm'))


@raises(ValueError)
@with_setup(setup_default_grid)
def test_name_setter_too_few_units():
    rmg.axis_units = ('m', )


@raises(ValueError)
@with_setup(setup_default_grid)
def test_name_setter_too_many_units():
    rmg.axis_units = ('m', 'cm', 'mm')

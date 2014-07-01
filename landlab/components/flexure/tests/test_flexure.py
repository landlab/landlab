#! /usr/bin/env python
"""
Unit tests for landlab.components.flexure.flexure
"""
from nose.tools import assert_equal, assert_true, assert_raises, with_setup
try:
    from nose.tools import assert_is_instance
except ImportError:
    from landlab.testing.tools import assert_is_instance
import numpy as np

from landlab import RasterModelGrid
from landlab.components.flexure.flexure import FlexureComponent


(_SHAPE, _SPACING, _ORIGIN) = ((20, 20), (10e3, 10e3), (0., 0.))
_ARGS = (_SHAPE, _SPACING, _ORIGIN)


def setup_grid():
    from landlab import RasterModelGrid
    grid = RasterModelGrid(20, 20, 10e3)
    flex = FlexureComponent(grid)
    globals().update({
        'flex': FlexureComponent(grid)
    })


@with_setup(setup_grid)
def test_name():
    assert_equal(flex.name, 'Flexure')


@with_setup(setup_grid)
def test_input_var_names():
    assert_equal(sorted(flex.input_var_names),
                 ['lithosphere__elevation',
                  'lithosphere__overlying_pressure',
                  'planet_surface_sediment__deposition_increment'])


@with_setup(setup_grid)
def test_output_var_names():
    assert_equal(sorted(flex.output_var_names),
                 ['lithosphere__elevation',
                  'lithosphere__elevation_increment'])


@with_setup(setup_grid)
def test_var_units():
    assert_equal(set(flex.input_var_names) |
                 set(flex.output_var_names),
                 set(flex.units))

    assert_equal(flex.units['lithosphere__elevation'], 'm')
    assert_equal(flex.units['lithosphere__elevation_increment'], 'm')
    assert_equal(flex.units['lithosphere__overlying_pressure'], 'Pa')
    assert_equal(
        flex.units['planet_surface_sediment__deposition_increment'], 'm')


@with_setup(setup_grid)
def test_grid_shape():
    assert_equal(flex.grid.number_of_node_rows, _SHAPE[0])
    assert_equal(flex.grid.number_of_node_columns, _SHAPE[1])


@with_setup(setup_grid)
def test_grid_xdimension():
    assert_equal(flex.grid.get_grid_xdimension(),
                 (_SHAPE[1] - 1) * _SPACING[1])


@with_setup(setup_grid)
def test_grid_ydimension():
    assert_equal(flex.grid.get_grid_xdimension(),
                 (_SHAPE[0] - 1) * _SPACING[0])


@with_setup(setup_grid)
def test_field_getters():
    for name in flex.grid['node']:
        field = flex.grid['node'][name]
        assert_is_instance(field, np.ndarray)
        assert_equal(field.shape,
                     (flex.grid.number_of_node_rows *
                      flex.grid.number_of_node_columns, ))

    assert_raises(KeyError, lambda: flex.grid['not_a_var_name'])


@with_setup(setup_grid)
def test_field_initialized_to_zero():
    for name in flex.grid['node']:
        field = flex.grid['node'][name]
        assert_true(np.all(field == 0.))

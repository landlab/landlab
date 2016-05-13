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
from landlab.components.flexure.flexure import Flexure


(_SHAPE, _SPACING, _ORIGIN) = ((20, 20), (10e3, 10e3), (0., 0.))
_ARGS = (_SHAPE, _SPACING, _ORIGIN)


def setup_grid():
    from landlab import RasterModelGrid
    grid = RasterModelGrid((20, 20), spacing=10e3)
    flex = Flexure(grid)
    globals().update({
        'flex': Flexure(grid)
    })


@with_setup(setup_grid)
def test_name():
    assert_equal(flex.name, 'Flexure')


@with_setup(setup_grid)
def test_input_var_names():
    assert_equal(flex.input_var_names,
                 ('lithosphere__overlying_pressure_increment', ))


@with_setup(setup_grid)
def test_output_var_names():
    assert_equal(flex.output_var_names,
                 ('lithosphere_surface__elevation_increment',))


@with_setup(setup_grid)
def test_var_units():
    assert_equal(set(flex.input_var_names) |
                 set(flex.output_var_names),
                 set(dict(flex.units).keys()))

    assert_equal(flex.var_units('lithosphere_surface__elevation_increment'), 'm')
    assert_equal(flex.var_units('lithosphere__overlying_pressure_increment'),
                 'Pa')


@with_setup(setup_grid)
def test_grid_shape():
    assert_equal(flex.grid.number_of_node_rows, _SHAPE[0])
    assert_equal(flex.grid.number_of_node_columns, _SHAPE[1])


@with_setup(setup_grid)
def test_grid_x_extent():
    assert_equal(flex.grid.extent[1], (_SHAPE[1] - 1) * _SPACING[1])


@with_setup(setup_grid)
def test_grid_y_extent():
    assert_equal(flex.grid.extent[0], (_SHAPE[0] - 1) * _SPACING[0])


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

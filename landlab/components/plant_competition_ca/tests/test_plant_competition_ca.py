"""
Unit tests for landlab.components.plant_competition_ca.plant_competition_ca
"""
from nose.tools import assert_equal, assert_true, assert_raises, with_setup
from numpy.testing import assert_array_almost_equal
try:
    from nose.tools import assert_is_instance
except ImportError:
    from landlab.testing.tools import assert_is_instance
import numpy as np

from landlab import RasterModelGrid
from landlab.components.plant_competition_ca.plant_competition_ca \
             import VegCA


(_SHAPE, _SPACING, _ORIGIN) = ((20, 20), (10e0, 10e0), (0., 0.))
_ARGS = (_SHAPE, _SPACING, _ORIGIN)


def setup_grid():
    from landlab import RasterModelGrid
    grid = RasterModelGrid((20, 20), spacing=10e0)
    ca_veg = VegCA(grid)
    globals().update({
        'ca_veg': VegCA(grid)
    })


@with_setup(setup_grid)
def test_name():
    assert_equal(ca_veg.name, 'Cellular Automata Plant Competition')


@with_setup(setup_grid)
def test_input_var_names():
    assert_equal(sorted(ca_veg.input_var_names),
                 ['vegetation__cumulative_water_stress',
                  'vegetation__plant_functional_type'])


@with_setup(setup_grid)
def test_output_var_names():
    assert_equal(sorted(ca_veg.output_var_names),
                 ['plant__age', 'plant__live_index'])


@with_setup(setup_grid)
def test_var_units():
    assert_equal(set(ca_veg.input_var_names) |
                 set(ca_veg.output_var_names),
                 set(dict(ca_veg.units).keys()))

    assert_equal(ca_veg.var_units('vegetation__cumulative_water_stress'),
                 'None')
    assert_equal(ca_veg.var_units('vegetation__plant_functional_type'),
                 'None')
    assert_equal(ca_veg.var_units('plant__live_index'), 'None')
    assert_equal(ca_veg.var_units('plant__age'), 'Years')


@with_setup(setup_grid)
def test_grid_shape():
    assert_equal(ca_veg.grid.number_of_node_rows, _SHAPE[0])
    assert_equal(ca_veg.grid.number_of_node_columns, _SHAPE[1])


@with_setup(setup_grid)
def test_grid_x_extent():
    assert_equal(ca_veg.grid.extent[1], (_SHAPE[1] - 1) * _SPACING[1])


@with_setup(setup_grid)
def test_grid_y_extent():
    assert_equal(ca_veg.grid.extent[0], (_SHAPE[0] - 1) * _SPACING[0])


@with_setup(setup_grid)
def test_field_getters():
    for name in ca_veg.grid['node']:
        field = ca_veg.grid['node'][name]
        assert_is_instance(field, np.ndarray)
        assert_equal(field.shape,
                     (ca_veg.grid.number_of_node_rows *
                      ca_veg.grid.number_of_node_columns, ))
                      
    for name in ca_veg.grid['cell']:
        field = ca_veg.grid['cell'][name]
        assert_is_instance(field, np.ndarray)
        assert_equal(field.shape,
                     (ca_veg.grid.number_of_cell_rows *
                      ca_veg.grid.number_of_cell_columns, ))

    assert_raises(KeyError, lambda: ca_veg.grid['not_a_var_name'])


@with_setup(setup_grid)
def test_field_initialized_to_zero():
    for name in ca_veg.grid['node']:
        field = ca_veg.grid['node'][name]
        assert_array_almost_equal(field, np.zeros(ca_veg.grid.number_of_nodes))

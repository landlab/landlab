"""
Unit tests for \
    landlab.components.vegetation_dynamics.vegetation_dynamics
"""
from nose.tools import assert_equal, assert_true, assert_raises, with_setup
from numpy.testing import assert_array_almost_equal
try:
    from nose.tools import assert_is_instance
except ImportError:
    from landlab.testing.tools import assert_is_instance
import numpy as np

from landlab import RasterModelGrid
from landlab.components.vegetation_dynamics.vegetation_dynamics \
             import Vegetation


(_SHAPE, _SPACING, _ORIGIN) = ((20, 20), (10e0, 10e0), (0., 0.))
_ARGS = (_SHAPE, _SPACING, _ORIGIN)


def setup_grid():
    from landlab import RasterModelGrid
    grid = RasterModelGrid((20, 20), spacing=10e0)
    grid['cell']['vegetation__plant_functional_type']= \
        np.zeros(grid.number_of_cells, dtype=int)
    Veg = Vegetation(grid)
    globals().update({
        'Veg': Vegetation(grid)
    })


@with_setup(setup_grid)
def test_name():
    assert_equal(Veg.name, 'Vegetation')


@with_setup(setup_grid)
def test_input_var_names():
    assert_equal(sorted(Veg.input_var_names),
                ['surface__evapotranspiration',
                 'surface__potential_evapotranspiration_30day_mean',
                 'surface__potential_evapotranspiration_rate',
                 'vegetation__plant_functional_type',
                 'vegetation__water_stress'])


@with_setup(setup_grid)
def test_output_var_names():
    assert_equal(sorted(Veg.output_var_names),
                ['vegetation__cover_fraction',
                 'vegetation__dead_biomass',
                 'vegetation__dead_leaf_area_index',
                 'vegetation__live_biomass',
                 'vegetation__live_leaf_area_index'])


@with_setup(setup_grid)
def test_var_units():
    assert_equal(set(Veg.input_var_names) |
                 set(Veg.output_var_names),
                 set(dict(Veg.units).keys()))

    assert_equal(Veg.var_units('vegetation__live_leaf_area_index'), 'None')
    assert_equal(Veg.var_units('vegetation__dead_leaf_area_index'), 'None')
    assert_equal(Veg.var_units('vegetation__cover_fraction'), 'None')
    assert_equal(Veg.var_units('surface__evapotranspiration'), 'mm')
    assert_equal(Veg.var_units('surface__potential_evapotranspiration_rate'),
                 'mm')
    assert_equal(Veg.var_units(
                'surface__potential_evapotranspiration_30day_mean'), 'mm')                 
    assert_equal(Veg.var_units('vegetation__water_stress'), 'None')
    assert_equal(Veg.var_units('vegetation__live_biomass'), 'g m^-2 d^-1')
    assert_equal(Veg.var_units('vegetation__dead_biomass'), 'g m^-2 d^-1')
    assert_equal(Veg.var_units('vegetation__plant_functional_type'), 'None')    


@with_setup(setup_grid)
def test_grid_shape():
    assert_equal(Veg.grid.number_of_node_rows, _SHAPE[0])
    assert_equal(Veg.grid.number_of_node_columns, _SHAPE[1])


@with_setup(setup_grid)
def test_grid_x_extent():
    assert_equal(Veg.grid.extent[1], (_SHAPE[1] - 1) * _SPACING[1])


@with_setup(setup_grid)
def test_grid_y_extent():
    assert_equal(Veg.grid.extent[0], (_SHAPE[0] - 1) * _SPACING[0])


@with_setup(setup_grid)
def test_field_getters():
    for name in Veg.grid['node']:
        field = Veg.grid['node'][name]
        assert_is_instance(field, np.ndarray)
        assert_equal(field.shape,
                     (Veg.grid.number_of_node_rows *
                      Veg.grid.number_of_node_columns, ))
                      
    for name in Veg.grid['cell']:
        field = Veg.grid['cell'][name]
        assert_is_instance(field, np.ndarray)
        assert_equal(field.shape,
                     (Veg.grid.number_of_cell_rows *
                      Veg.grid.number_of_cell_columns, ))

    assert_raises(KeyError, lambda: Veg.grid['not_a_var_name'])


@with_setup(setup_grid)
def test_field_initialized_to_zero():
    for name in Veg.grid['node']:
        field = Veg.grid['node'][name]
        assert_array_almost_equal(field, np.zeros(Veg.grid.number_of_nodes))
    for name in Veg.grid['cell']:
        field = Veg.grid['cell'][name]
        assert_array_almost_equal(field, np.zeros(Veg.grid.number_of_cells))
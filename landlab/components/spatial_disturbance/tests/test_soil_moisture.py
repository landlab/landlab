"""
Unit tests for \
    landlab.components.soil_moisture.soil_moisture_dynamics
"""
from nose.tools import assert_equal, assert_true, assert_raises, with_setup
from numpy.testing import assert_array_almost_equal
try:
    from nose.tools import assert_is_instance
except ImportError:
    from landlab.testing.tools import assert_is_instance
import numpy as np

from landlab import RasterModelGrid
from landlab.components.soil_moisture.soil_moisture_dynamics \
             import SoilMoisture


(_SHAPE, _SPACING, _ORIGIN) = ((20, 20), (10e0, 10e0), (0., 0.))
_ARGS = (_SHAPE, _SPACING, _ORIGIN)


def setup_grid():
    from landlab import RasterModelGrid
    grid = RasterModelGrid((20, 20), spacing=10e0)
    grid['cell']['vegetation__plant_functional_type']= \
        np.zeros(grid.number_of_cells, dtype=int)
    SM = SoilMoisture(grid)
    globals().update({
        'SM': SoilMoisture(grid)
    })


@with_setup(setup_grid)
def test_name():
    assert_equal(SM.name, 'Soil Moisture')


@with_setup(setup_grid)
def test_input_var_names():
    assert_equal(sorted(SM.input_var_names),
                ['rainfall__daily_depth',
                 'soil_moisture__initial_saturation_fraction',
                 'surface__potential_evapotranspiration_rate',
                 'vegetation__cover_fraction',
                 'vegetation__live_leaf_area_index',
                 'vegetation__plant_functional_type'])


@with_setup(setup_grid)
def test_output_var_names():
    assert_equal(sorted(SM.output_var_names),
                ['soil_moisture__root_zone_leakage',
                 'soil_moisture__saturation_fraction',
                 'surface__evapotranspiration',
                 'surface__runoff',
                 'vegetation__water_stress',])


@with_setup(setup_grid)
def test_var_units():
    assert_equal(set(SM.input_var_names) |
                 set(SM.output_var_names),
                 set(dict(SM.units).keys()))

    assert_equal(SM.var_units('vegetation__cover_fraction'), 'None')
    assert_equal(SM.var_units('vegetation__live_leaf_area_index'), 'None')
    assert_equal(SM.var_units('surface__potential_evapotranspiration_rate'),
                 'mm')
    assert_equal(SM.var_units('vegetation__plant_functional_type'), 'None')
    assert_equal(SM.var_units('vegetation__water_stress'), 'None')
    assert_equal(SM.var_units('soil_moisture__saturation_fraction'), 'None')
    assert_equal(SM.var_units('soil_moisture__initial_saturation_fraction'),
                 'None')
    assert_equal(SM.var_units('soil_moisture__root_zone_leakage'), 'mm')
    assert_equal(SM.var_units('surface__runoff'), 'mm')
    assert_equal(SM.var_units('surface__evapotranspiration'), 'mm')
    assert_equal(SM.var_units('rainfall__daily_depth'), 'mm')


@with_setup(setup_grid)
def test_grid_shape():
    assert_equal(SM.grid.number_of_node_rows, _SHAPE[0])
    assert_equal(SM.grid.number_of_node_columns, _SHAPE[1])


@with_setup(setup_grid)
def test_grid_x_extent():
    assert_equal(SM.grid.extent[1], (_SHAPE[1] - 1) * _SPACING[1])


@with_setup(setup_grid)
def test_grid_y_extent():
    assert_equal(SM.grid.extent[0], (_SHAPE[0] - 1) * _SPACING[0])


@with_setup(setup_grid)
def test_field_getters():
    for name in SM.grid['node']:
        field = SM.grid['node'][name]
        assert_is_instance(field, np.ndarray)
        assert_equal(field.shape,
                     (SM.grid.number_of_node_rows *
                      SM.grid.number_of_node_columns, ))
                      
    for name in SM.grid['cell']:
        field = SM.grid['cell'][name]
        assert_is_instance(field, np.ndarray)
        assert_equal(field.shape,
                     (SM.grid.number_of_cell_rows *
                      SM.grid.number_of_cell_columns, ))

    assert_raises(KeyError, lambda: SM.grid['not_a_var_name'])


@with_setup(setup_grid)
def test_field_initialized_to_zero():
    for name in SM.grid['node']:
        field = SM.grid['node'][name]
        assert_array_almost_equal(field, np.zeros(SM.grid.number_of_nodes))
    for name in SM.grid['cell']:
        field = SM.grid['cell'][name]
        assert_array_almost_equal(field, np.zeros(SM.grid.number_of_cells))
"""
Unit tests for landlab.components.pet.potential_evapotranspiration_field
"""
import pytest
from nose.tools import with_setup
from numpy.testing import assert_array_almost_equal
import numpy as np

from landlab import RasterModelGrid
from landlab.components.pet.potential_evapotranspiration_field \
             import PotentialEvapotranspiration


(_SHAPE, _SPACING, _ORIGIN) = ((20, 20), (10e0, 10e0), (0., 0.))
_ARGS = (_SHAPE, _SPACING, _ORIGIN)


def setup_grid():
    from landlab import RasterModelGrid
    grid = RasterModelGrid((20, 20), spacing=10e0)
    PET = PotentialEvapotranspiration(grid)
    globals().update({
        'PET': PotentialEvapotranspiration(grid)
    })


@with_setup(setup_grid)
def test_name():
    assert PET.name == 'Potential Evapotranspiration'


@with_setup(setup_grid)
def test_input_var_names():
    assert PET.input_var_names == ('radiation__ratio_to_flat_surface', )


@with_setup(setup_grid)
def test_output_var_names():
    assert sorted(PET.output_var_names) == [
        'radiation__incoming_shortwave_flux',
        'radiation__net_flux',
        'radiation__net_longwave_flux',
        'radiation__net_shortwave_flux',
        'surface__potential_evapotranspiration_rate',
    ]


@with_setup(setup_grid)
def test_var_units():
    assert set(PET.input_var_names) | set(PET.output_var_names) == set(dict(PET.units).keys())

    assert PET.var_units('radiation__incoming_shortwave_flux') == 'W/m^2'
    assert PET.var_units('radiation__net_flux') == 'W/m^2'
    assert PET.var_units('radiation__net_longwave_flux') == 'W/m^2'
    assert PET.var_units('radiation__net_shortwave_flux') == 'W/m^2'
    assert PET.var_units('radiation__ratio_to_flat_surface') == 'None'
    assert PET.var_units('surface__potential_evapotranspiration_rate') == 'mm'


@with_setup(setup_grid)
def test_grid_shape():
    assert PET.grid.number_of_node_rows == _SHAPE[0]
    assert PET.grid.number_of_node_columns == _SHAPE[1]


@with_setup(setup_grid)
def test_grid_x_extent():
    assert PET.grid.extent[1] == (_SHAPE[1] - 1) * _SPACING[1]


@with_setup(setup_grid)
def test_grid_y_extent():
    assert PET.grid.extent[0] == (_SHAPE[0] - 1) * _SPACING[0]


@with_setup(setup_grid)
def test_field_getters():
    for name in PET.grid['node']:
        field = PET.grid['node'][name]
        assert isinstance(field, np.ndarray)
        assert field.shape == (PET.grid.number_of_node_rows *
                               PET.grid.number_of_node_columns, )
                      
    for name in PET.grid['cell']:
        field = PET.grid['cell'][name]
        assert isinstance(field, np.ndarray)
        assert field.shape == (PET.grid.number_of_cell_rows *
                               PET.grid.number_of_cell_columns, )

    with pytest.raises(KeyError):
        PET.grid['not_a_var_name']


@with_setup(setup_grid)
def test_field_initialized_to_zero():
    for name in PET.grid['node']:
        field = PET.grid['node'][name]
        assert_array_almost_equal(field, np.zeros(PET.grid.number_of_nodes))
    for name in PET.grid['cell']:
        field = PET.grid['cell'][name]
        assert_array_almost_equal(field, np.zeros(PET.grid.number_of_cells))

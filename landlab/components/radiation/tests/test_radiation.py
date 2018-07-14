"""
Unit tests for landlab.components.radiation.radiation
"""
import pytest
from nose.tools import with_setup
from numpy.testing import assert_array_almost_equal
import numpy as np

from landlab import RasterModelGrid
from landlab.components.radiation.radiation import Radiation


(_SHAPE, _SPACING, _ORIGIN) = ((20, 20), (10e0, 10e0), (0., 0.))
_ARGS = (_SHAPE, _SPACING, _ORIGIN)


def setup_grid():
    from landlab import RasterModelGrid
    grid = RasterModelGrid((20, 20), spacing=10e0)
    rad = Radiation(grid)
    globals().update({
        'rad': Radiation(grid)
    })


@with_setup(setup_grid)
def test_name():
    assert rad.name == 'Radiation'


@with_setup(setup_grid)
def test_input_var_names():
    assert rad.input_var_names == ('topographic__elevation',)


@with_setup(setup_grid)
def test_output_var_names():
    assert sorted(rad.output_var_names) == [
        'radiation__incoming_shortwave_flux',
        'radiation__net_shortwave_flux',
        'radiation__ratio_to_flat_surface'
    ]


@with_setup(setup_grid)
def test_var_units():
    assert set(rad.input_var_names) | set(rad.output_var_names) == set(dict(rad.units).keys())

    assert rad.var_units('topographic__elevation') == 'm'
    assert rad.var_units('radiation__incoming_shortwave_flux') == 'W/m^2'
    assert rad.var_units('radiation__net_shortwave_flux') == 'W/m^2'
    assert rad.var_units('radiation__ratio_to_flat_surface') == 'None'


@with_setup(setup_grid)
def test_grid_shape():
    assert rad.grid.number_of_node_rows == _SHAPE[0]
    assert rad.grid.number_of_node_columns == _SHAPE[1]


@with_setup(setup_grid)
def test_grid_x_extent():
    assert rad.grid.extent[1] == (_SHAPE[1] - 1) * _SPACING[1]


@with_setup(setup_grid)
def test_grid_y_extent():
    assert rad.grid.extent[0] == (_SHAPE[0] - 1) * _SPACING[0]


@with_setup(setup_grid)
def test_field_getters():
    for name in rad.grid['node']:
        field = rad.grid['node'][name]
        assert isinstance(field, np.ndarray)
        assert field.shape == (rad.grid.number_of_node_rows * rad.grid.number_of_node_columns, )
                      
    for name in rad.grid['cell']:
        field = rad.grid['cell'][name]
        assert isinstance(field, np.ndarray)
        assert field.shape == (rad.grid.number_of_cell_rows * rad.grid.number_of_cell_columns, )

    with pytest.raises(KeyError):
        rad.grid['not_a_var_name']


@with_setup(setup_grid)
def test_field_initialized_to_zero():
    for name in rad.grid['node']:
        field = rad.grid['node'][name]
        assert_array_almost_equal(field, np.zeros(rad.grid.number_of_nodes))
    for name in rad.grid['cell']:
        if name == 'Slope' or name == 'Aspect':
            continue
        field = rad.grid['cell'][name]
        assert_array_almost_equal(field, np.zeros(rad.grid.number_of_cells))

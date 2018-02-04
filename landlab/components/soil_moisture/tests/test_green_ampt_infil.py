# -*- coding: utf-8 -*-
"""
Unit tests for landlab.components.soil_moisture.SoilInfiltrationGreenAmpt

last updated: 3/14/16
"""
from nose.tools import assert_equal, assert_true, assert_raises, with_setup
try:
    from nose.tools import assert_is_instance
except ImportError:
    from landlab.testing.tools import assert_is_instance
import numpy as np

from landlab import RasterModelGrid
from landlab.components.soil_moisture import SoilInfiltrationGreenAmpt

(_SHAPE, _SPACING, _ORIGIN) = ((10, 10), (25, 25), (0., 0.))
_ARGS = (_SHAPE, _SPACING, _ORIGIN)


def setup_grid():
    from landlab import RasterModelGrid
    grid = RasterModelGrid((10, 10), spacing=25)
    grid.add_ones('node', 'soil_water_infiltration__depth', dtype=float)
    grid.add_ones('node', 'surface_water__depth')
    hydraulic_conductivity = (2.5 *  (10**-5))
    grid['node']['surface_water__depth'] *= 0.5
    grid['node']['soil_water_infiltration__depth'] *= (10**-5)
    SI = SoilInfiltrationGreenAmpt(grid,
        hydraulic_conductivity=hydraulic_conductivity,
         soil_bulk_density=1700., rock_density=2650.,
         initial_soil_moisture_content=0.2, soil_type='silt loam',
         volume_fraction_coarse_fragments=0.6,
         coarse_sed_flag=False,
         surface_water_minimum_depth=1.e-7,
         soil_pore_size_distribution_index=None,
         soil_bubbling_pressure=None,
         wetting_front_capillary_pressure_head=None)

    globals().update({
        'SI': SoilInfiltrationGreenAmpt(grid)})


@with_setup(setup_grid)
def test_SI_name():
    assert_equal(SI.name, 'SoilInfiltrationGreenAmpt')


@with_setup(setup_grid)
def test_SI_input_var_names():
    assert_equal(SI.input_var_names,  ('surface_water__depth',
                                          'soil_water_infiltration__depth',))


@with_setup(setup_grid)
def test_SI_output_var_names():
    assert_equal(SI.output_var_names, ('surface_water__depth',
                                          'soil_water_infiltration__depth',))

@with_setup(setup_grid)
def test_SI_var_units():
    assert_equal(set(SI.input_var_names) |
                 set(SI.output_var_names),
                 set(dict(SI.units).keys()))

    assert_equal(SI.var_units('surface_water__depth'), 'm')
    assert_equal(SI.var_units('soil_water_infiltration__depth'), 'm')


@with_setup(setup_grid)
def test_grid_shape():
    assert_equal(SI.grid.number_of_node_rows, _SHAPE[0])
    assert_equal(SI.grid.number_of_node_columns, _SHAPE[1])


@with_setup(setup_grid)
def test_calc_soil_pressure():
    np.testing.assert_almost_equal(SI.calc_soil_pressure('silt loam'),
                            0.1647870740305523, decimal=6)
@with_setup(setup_grid)
def test_calc_soil_head():
    soil_props = SoilInfiltrationGreenAmpt.SOIL_PROPS['loam']
    np.testing.assert_almost_equal(SI.calc_pressure_head(soil_props[0],
            soil_props[1]), 0.087498292, decimal=6)

@with_setup(setup_grid)
def test_calc_moisture_deficit():
    np.testing.assert_almost_equal(SI.calc_moisture_deficit(
        soil_bulk_density=1700., rock_density=2650.,
        volume_fraction_coarse_fragments=0.,
        soil_moisture_content=0.2), 0.15849056603, decimal=6)

def test_run_one_step():
    from landlab import RasterModelGrid
    grid = RasterModelGrid((10, 10), spacing=25)
    grid.add_ones('node', 'soil_water_infiltration__depth', dtype=float)
    grid.add_ones('node', 'surface_water__depth')
    hydraulic_conductivity = (2.5 *  (10**-6))
    grid['node']['surface_water__depth'] *= 5.0
    grid['node']['soil_water_infiltration__depth'] *= (10**-5)
    SI = SoilInfiltrationGreenAmpt(grid,
        hydraulic_conductivity=hydraulic_conductivity,
         soil_bulk_density=1700., rock_density=2650.,
         initial_soil_moisture_content=0.2, soil_type='silt loam',
         volume_fraction_coarse_fragments=0.6,
         coarse_sed_flag=False,
         surface_water_minimum_depth=1.e-7,
         soil_pore_size_distribution_index=None,
         soil_bubbling_pressure=None,
         wetting_front_capillary_pressure_head=None)

    SI.run_one_step(dt=5)
    np.testing.assert_almost_equal(grid['node']['surface_water__depth'][0],
                3.97677483519, decimal=6)

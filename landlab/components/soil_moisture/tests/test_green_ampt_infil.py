# -*- coding: utf-8 -*-
"""
Unit tests for landlab.components.soil_moisture.SoilInfiltrationGreenAmpt

last updated: 3/14/16
"""
import numpy as np

from landlab import RasterModelGrid
from landlab.components.soil_moisture import SoilInfiltrationGreenAmpt

(_SHAPE, _SPACING, _ORIGIN) = ((10, 10), (25, 25), (0.0, 0.0))
_ARGS = (_SHAPE, _SPACING, _ORIGIN)


def test_SI_name(si):
    assert si.name == "SoilInfiltrationGreenAmpt"


def test_SI_input_var_names(si):
    assert si.input_var_names == (
        "surface_water__depth",
        "soil_water_infiltration__depth",
    )


def test_SI_output_var_names(si):
    assert si.output_var_names == (
        "surface_water__depth",
        "soil_water_infiltration__depth",
    )


def test_SI_var_units(si):
    assert set(si.input_var_names) | set(si.output_var_names) == set(
        dict(si.units).keys()
    )

    assert si.var_units("surface_water__depth") == "m"
    assert si.var_units("soil_water_infiltration__depth") == "m"


def test_grid_shape(si):
    assert si.grid.number_of_node_rows == _SHAPE[0]
    assert si.grid.number_of_node_columns == _SHAPE[1]


def test_calc_soil_pressure(si):
    np.testing.assert_almost_equal(
        si.calc_soil_pressure("silt loam"), 0.1647870740305523, decimal=6
    )


def test_calc_soil_head(si):
    soil_props = SoilInfiltrationGreenAmpt.SOIL_PROPS["loam"]
    np.testing.assert_almost_equal(
        si.calc_pressure_head(soil_props[0], soil_props[1]), 0.087498292, decimal=6
    )


def test_calc_moisture_deficit(si):
    np.testing.assert_almost_equal(
        si.calc_moisture_deficit(
            soil_bulk_density=1700.0,
            rock_density=2650.0,
            volume_fraction_coarse_fragments=0.0,
            soil_moisture_content=0.2,
        ),
        0.15849056603,
        decimal=6,
    )


def test_run_one_step():
    grid = RasterModelGrid((10, 10), xy_spacing=25)
    grid.add_ones("node", "soil_water_infiltration__depth", dtype=float)
    grid.add_ones("node", "surface_water__depth")
    hydraulic_conductivity = 2.5 * (10 ** -6)
    grid["node"]["surface_water__depth"] *= 5.0
    grid["node"]["soil_water_infiltration__depth"] *= 10 ** -5
    SI = SoilInfiltrationGreenAmpt(
        grid,
        hydraulic_conductivity=hydraulic_conductivity,
        soil_bulk_density=1700.0,
        rock_density=2650.0,
        initial_soil_moisture_content=0.2,
        soil_type="silt loam",
        volume_fraction_coarse_fragments=0.6,
        coarse_sed_flag=False,
        surface_water_minimum_depth=1.0e-7,
        soil_pore_size_distribution_index=None,
        soil_bubbling_pressure=None,
        wetting_front_capillary_pressure_head=None,
    )

    SI.run_one_step(dt=5)
    np.testing.assert_almost_equal(
        grid["node"]["surface_water__depth"][0], 3.97677483519, decimal=6
    )

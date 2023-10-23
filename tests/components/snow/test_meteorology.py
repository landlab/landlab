"""
Unit tests for landlab.components.snow.meteorology

@author Tian Gan  Sept 2023
"""

import pytest
from datetime import datetime
import numpy as np
from numpy.testing import assert_almost_equal

from landlab import RasterModelGrid
from landlab.components import Meteorology


def test_create_fields():
    """check to make sure the right fields are created"""

    grid = RasterModelGrid((2, 2))
    grid.add_full("atmosphere_bottom_air__temperature", 1, at="node")
    grid.add_full("land_surface__temperature", -1, at="node")
    grid.add_full("land_surface__latitude", 40, at="node")
    grid.add_full("land_surface__longitude", -105, at="node")

    Meteorology(grid, start_datetime='2023-01-01 12:00:00')

    assert len(grid.at_node.keys()) == 35
    assert_almost_equal(grid.at_node["atmosphere_bottom_air__temperature"], 1)
    assert_almost_equal(grid.at_node["land_surface__temperature"], -1)
    assert_almost_equal(grid.at_node["land_surface__latitude"], 40)
    assert_almost_equal(grid.at_node["land_surface__longitude"], -105)
    assert_almost_equal(grid.at_node["land_surface__aspect_angle"], 0)
    assert_almost_equal(grid.at_node["land_surface__slope_angle"], 0)
    assert_almost_equal(grid.at_node["snowpack__depth"], 0)
    assert_almost_equal(grid.at_node["land_surface__albedo"], 0.3)
    assert_almost_equal(grid.at_node["land_surface__emissivity"], 0.98)
    assert_almost_equal(
        grid.at_node["atmosphere_aerosol_dust__reduction_of_transmittance"], 0.08
    )
    assert_almost_equal(
        grid.at_node["atmosphere_bottom_air__brutsaert_emissivity_canopy_factor"], 0
    )
    assert_almost_equal(
        grid.at_node["atmosphere_bottom_air__brutsaert_emissivity_cloud_factor"], 0
    )
    assert_almost_equal(
        grid.at_node["atmosphere_bottom_air_water-vapor__relative_saturation"], 0.5
    )
    assert_almost_equal(
        grid.at_node["atmosphere_bottom_air__pressure"], 1000
    )
    assert_almost_equal(
        grid.at_node["atmosphere_bottom_air_flow__log_law_roughness_length"], 0.02
    )
    assert_almost_equal(
        grid.at_node["atmosphere_bottom_air_flow__speed_reference_height"], 10
    )
    assert_almost_equal(
        grid.at_node["atmosphere_bottom_air_flow__reference-height_speed"], 3
    )


def test_assign_parameters():
    """Test when parameters are updated"""

    grid = RasterModelGrid((2, 2))
    grid.add_full("atmosphere_bottom_air__temperature", 1, at="node")
    grid.add_full("land_surface__temperature", -1, at="node")
    grid.add_full("land_surface__latitude", 40, at="node")
    grid.add_full("land_surface__longitude", -105, at="node")

    with pytest.raises(ValueError):
        met = Meteorology(grid, start_datetime='23/01/0112:00:00')

    met = Meteorology(grid, start_datetime='2023-01-01 12:00:00')

    # constants
    assert met._g == np.float64(9.81)
    assert met._kappa == np.float64(0.408)
    assert met._Lv == np.float64(2500000)
    assert met._Lf == np.float64(334000)
    assert met._sigma == np.float64(5.67e-8)
    assert met._C_to_K == np.float64(273.15)

    assert met._one_seventh == np.float64(1) / 7
    assert met._hours_per_day == np.float64(24)
    assert met._latent_heat_constant == np.float64(0.662)

    # parameters
    assert met._datetime_obj == datetime(2023, 1, 1, 12, 0)
    assert met._GMT_offset == 0
    assert met._rho_H2O == 1000
    assert met._rho_air == 1.2614
    assert met._Cp_air == 1005.7
    assert not met._satterlund

    with pytest.raises(AssertionError):
        met.rho_H2O = -1

    with pytest.raises(AssertionError):
        met.rho_air = -1

    with pytest.raises(AssertionError):
        met.Cp_air = -1


# def test_snow_accumulation():
#     """Test when there is only snow accumulation"""
#
#     grid = RasterModelGrid((2, 2))
#     grid.add_full("atmosphere_water__precipitation_leq-volume_flux", 0.001, at="node")
#     grid.add_full("atmosphere_bottom_air__temperature", -1, at="node")
#     grid.add_full("land_surface__temperature", -1, at="node")
#     grid.add_full("land_surface_net-total-energy__energy_flux", 0, at="node")
#     grid.add_full("snowpack__liquid-equivalent_depth", 1, at="node")
#     grid.add_full("snowpack__z_mean_of_mass-per-volume_density", 200, at="node")
#     grid.add_full(
#         "snowpack__z_mean_of_mass-specific_isobaric_heat_capacity", 2000, at="node"
#     )
#
#     sm = SnowEnergyBalance(grid, grid_area=100)
#     assert sm.vol_SM == 0, f"wrong vol_SM value is {sm.vol_SM}"
#     assert sm.vol_swe == 400, f"wrong vol_swe value is {sm.vol_swe}"
#     assert_almost_equal(grid.at_node["snowpack__energy-per-area_cold_content"], 2e6)
#
#     sm.run_one_step(1000)  # dt = 1000 sec
#     assert_almost_equal(grid.at_node["snowpack__depth"], 10)
#     assert_almost_equal(grid.at_node["snowpack__melt_volume_flux"], 0)
#     # TODO: may need to change Ecc using new code for update_code_content
#     assert_almost_equal(grid.at_node["snowpack__energy-per-area_cold_content"], 2e6)
#     assert_almost_equal(
#         grid.at_node["atmosphere_water__time_integral_snowfall_leq-volume_flux"], 1
#     )
#     assert_almost_equal(grid.at_node["snowpack__time_integral_melt_volume_flux"], 0)
#
#     assert sm.vol_SM == 0, f"wrong vol_SM value is {sm.vol_SM}"
#     assert sm.vol_swe == 800, f"wrong vol_swe value is {sm.vol_swe}"





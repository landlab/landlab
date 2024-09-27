"""
Unit tests for landlab.components.snow.meteorology

@author Tian Gan  Sept 2023
"""

from datetime import datetime

import numpy as np
import pytest
from numpy.testing import assert_almost_equal

from landlab import RasterModelGrid
from landlab.components import Meteorology

PARAMS = {
    "test_case_1": {
        "air_temp": 8.0,
        "surf_temp": 4.0,
        "wind_speed": 2.0,
        "richardson_number": 0.34892406,
        "aero_conductance": 0.00192022,
        "dew_point": -2.5,
        "h_swe": 1,
        "air_pressure": 1000,
        "air_vapor_pressure": 5.384721,
        "surf_vapor_pressure": 2.971155,
        "canopy_factor": 0.98,
        "cloud_factor": 0.62,
        "air_emissivity": 0.99526025,
        "surf_emissivity": 0.98,
        "tsn_offset": 0.07271961,
        "latitude": 40,
        "alpha": 0.5,
        "beta": 0.2,
        "prec_water": 0.1,
        "albedo": 0.3,
        "dust_atten": 0.1,
    },
    "test_case_2": {
        "air_temp": [8.0],
        "surf_temp": [4.0],
        "wind_speed": [2.0],
        "richardson_number": [0.34892406],
        "aero_conductance": [0.00192022],
        "dew_point": [-2.5],
        "h_swe": [1],
        "air_pressure": [1000],
        "air_vapor_pressure": [5.384721],
        "surf_vapor_pressure": [2.971155],
        "canopy_factor": [0.98],
        "cloud_factor": [0.62],
        "air_emissivity": [0.99526025],
        "surf_emissivity": [0.98],
        "tsn_offset": [0.07271961],
        "latitude": [40],
        "alpha": [0.5],
        "beta": [0.2],
        "prec_water": [0.1],
        "albedo": [0.3],
        "dust_atten": [0.1],
    },
    "test_case_3": {
        "air_temp": [8.0, 8.0],
        "surf_temp": [4.0, 4.0],
        "wind_speed": [2.0, 2.0],
        "richardson_number": [0.348924, 0.348924],
        "aero_conductance": [0.00192022, 0.00192022],
        "dew_point": [-2.5, -2.5],
        "h_swe": [1, 1],
        "air_pressure": [1000, 1000],
        "air_vapor_pressure": [5.384721, 5.384721],
        "surf_vapor_pressure": [2.971155, 2.971155],
        "canopy_factor": [0.98, 0.98],
        "cloud_factor": [0.62, 0.62],
        "air_emissivity": [0.99526025, 0.99526025],
        "surf_emissivity": [0.98, 0.98],
        "tsn_offset": [0.07271961, 0.07271961],
        "latitude": [40, 40],
        "alpha": [0.5, 0.5],
        "beta": [0.2, 0.2],
        "prec_water": [0.1, 0.1],
        "albedo": [0.3, 0.3],
        "dust_atten": [0.1, 0.1],
    },
    "test_case_4": {
        "air_temp": [[8.0, 8.0], [8.0, 8.0]],
        "surf_temp": [[4.0, 4.0], [4.0, 4.0]],
        "wind_speed": [[2.0, 2.0], [2.0, 2.0]],
        "richardson_number": [[0.348924, 0.348924], [0.348924, 0.348924]],
        "aero_conductance": [[0.00192022, 0.00192022], [0.00192022, 0.00192022]],
        "dew_point": [[-2.5, -2.5], [-2.5, -2.5]],
        "h_swe": [[1, 1], [1, 1]],
        "air_pressure": [[1000, 1000], [1000, 1000]],
        "air_vapor_pressure": [[5.384721, 5.384721], [5.384721, 5.384721]],
        "surf_vapor_pressure": [[2.971155, 2.971155], [2.971155, 2.971155]],
        "canopy_factor": [[0.98, 0.98], [0.98, 0.98]],
        "cloud_factor": [[0.62, 0.62], [0.62, 0.62]],
        "air_emissivity": [[0.99526025, 0.99526025], [0.99526025, 0.99526025]],
        "surf_emissivity": [[0.98, 0.98], [0.98, 0.98]],
        "tsn_offset": [[0.07271961, 0.07271961], [0.07271961, 0.07271961]],
        "latitude": [[40, 40], [40, 40]],
        "alpha": [[0.5, 0.5], [0.5, 0.5]],
        "beta": [[0.2, 0.2], [0.2, 0.2]],
        "prec_water": [[0.1, 0.1], [0.1, 0.1]],
        "albedo": [[0.3, 0.3], [0.3, 0.3]],
        "dust_atten": [[0.1, 0.1], [0.1, 0.1]],
    },
}


def test_create_fields():
    """Test to create all fields"""

    grid = RasterModelGrid((2, 2))
    grid.add_full("atmosphere_bottom_air__temperature", 1, at="node")
    grid.add_full("land_surface__temperature", -1, at="node")
    grid.add_full("land_surface__latitude", 40, at="node")
    grid.add_full("land_surface__longitude", -105, at="node")

    Meteorology(grid, start_datetime="2023-01-01 12:00:00", calc_input=True)

    assert len(grid.at_node.keys()) == 30
    assert_almost_equal(grid.at_node["atmosphere_bottom_air__temperature"], 1)
    assert_almost_equal(grid.at_node["land_surface__latitude"], 40)
    assert_almost_equal(grid.at_node["land_surface__longitude"], -105)
    assert_almost_equal(grid.at_node["land_surface__aspect_angle"], 0)
    assert_almost_equal(grid.at_node["land_surface__slope_angle"], 0)
    assert_almost_equal(grid.at_node["snowpack__liquid-equivalent_depth"], 0)
    assert_almost_equal(grid.at_node["land_surface__albedo"], 0.3)
    assert_almost_equal(grid.at_node["land_surface__emissivity"], 0.98)
    assert_almost_equal(
        grid.at_node["atmosphere_aerosol_dust__reduction_of_transmittance"], 0.0
    )
    assert_almost_equal(
        grid.at_node["atmosphere_bottom_air__brutsaert_emissivity_canopy_factor"], 0.0
    )
    assert_almost_equal(
        grid.at_node["atmosphere_bottom_air__brutsaert_emissivity_cloud_factor"], 0.0
    )
    assert_almost_equal(
        grid.at_node["atmosphere_bottom_air_water-vapor__relative_saturation"], 0.5
    )
    assert_almost_equal(grid.at_node["atmosphere_bottom_air__pressure"], 1013.25)
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
        Meteorology(grid, start_datetime="2023/01/01")

    met = Meteorology(grid, start_datetime="2023-01-01 12:00:00")

    # parameters
    assert met._datetime_obj == datetime(2023, 1, 1, 12, 0)
    assert met.gmt_offset == 0
    assert met.rho_air == 1.2614
    assert met.cp_air == 1005.7
    assert met.roughness_length == 0.02
    assert met.reference_height == 10
    assert not met.satterlund
    assert not met.calc_input
    assert not met.clear_sky

    with pytest.raises(ValueError):
        met.gmt_offset = -13

    with pytest.raises(ValueError):
        met.rho_air = -1

    with pytest.raises(ValueError):
        met.cp_air = -1

    with pytest.raises(ValueError):
        met.roughness_length = -1

    with pytest.raises(ValueError):
        met.reference_height = -1

    with pytest.raises(ValueError):
        met.satterlund = -1

    with pytest.raises(ValueError):
        met.calc_input = -1

    with pytest.raises(ValueError):
        met.clear_sky = -1

@pytest.mark.parametrize(
    "params",
    [
        PARAMS["test_case_1"],
        PARAMS["test_case_2"],
        PARAMS["test_case_3"],
        PARAMS["test_case_4"],
    ],
)
def test_calc_bulk_richardson_number(params):
    """Test the calculation of the bulk Richardson number."""

    richardson_number = Meteorology.calc_bulk_richardson_number(
        reference_height=10,
        air_temp=params["air_temp"],
        surf_temp=params["surf_temp"],
        wind_speed=params["wind_speed"],
    )

    assert np.allclose(richardson_number, 0.3489240)


@pytest.mark.parametrize(
    "params",
    [
        PARAMS["test_case_1"],
        PARAMS["test_case_2"],
        PARAMS["test_case_3"],
        PARAMS["test_case_4"],
    ],
)
def test_calc_bulk_aero_conductance(params):
    """Test the calculation of the bulk aero conductance."""

    aero_conductance = Meteorology.calc_bulk_aero_conductance(
        reference_height=10,
        roughness_length=0.02,
        richardson_number=params["richardson_number"],
        wind_speed=params["wind_speed"],
        air_temp=params["air_temp"],
        surf_temp=params["surf_temp"],
    )
    assert np.allclose(aero_conductance, 0.00192022)


@pytest.mark.parametrize(
    "params",
    [
        PARAMS["test_case_1"],
        PARAMS["test_case_2"],
        PARAMS["test_case_3"],
        PARAMS["test_case_4"],
    ],
)
def test_calc_sensible_heat_flux(params):
    """test net sensible heat flux"""

    sensible_heat_flux = Meteorology.calc_sensible_heat_flux(
        rho_air=1.2614,
        cp_air=1005.7,
        aero_conductance=params["aero_conductance"],
        air_temp=params["air_temp"],
        surf_temp=params["surf_temp"],
    )

    assert np.allclose(sensible_heat_flux, 9.7438)


@pytest.mark.parametrize(
    "params",
    [
        PARAMS["test_case_1"],
        PARAMS["test_case_2"],
        PARAMS["test_case_3"],
        PARAMS["test_case_4"],
    ],
)
def test_calc_saturation_vapor_pressure(params):
    """Test to calculate saturation vapor pressure"""

    # Brutsaert method
    sat_vapor_pressure = Meteorology.calc_saturation_vapor_pressure(
        temp=params["air_temp"], satterlund=False, millibar=True
    )
    assert np.allclose(sat_vapor_pressure, 10.74170541)

    # Satterlund method
    sat_vapor_pressure = Meteorology.calc_saturation_vapor_pressure(
        temp=params["surf_temp"], satterlund=True, millibar=False
    )
    assert np.allclose(sat_vapor_pressure, 0.81285415)


@pytest.mark.parametrize(
    "sat_vapor_pressure", [10, [10.0], [10, 10], [[10, 10], [10, 10]]]
)
def test_calc_vapor_pressure(sat_vapor_pressure):
    """Test to calculate vapor pressure"""

    vapor_pressure = Meteorology.calc_vapor_pressure(
        sat_vapor_pressure=sat_vapor_pressure, relative_humidity=0.5
    )

    assert np.allclose(vapor_pressure, 5)


@pytest.mark.parametrize("air_vapor_pressure", [5, [5.0], [5, 5], [[5, 5], [5, 5]]])
def test_calc_dew_point(air_vapor_pressure):
    """Test to calculate dew point"""

    dew_point = Meteorology.calc_dew_point(air_vapor_pressure)

    assert np.allclose(dew_point, -2.73544824)


@pytest.mark.parametrize(
    "params",
    [
        PARAMS["test_case_1"],
        PARAMS["test_case_2"],
        PARAMS["test_case_3"],
        PARAMS["test_case_4"],
    ],
)
def test_calc_surf_temp(params):
    """Test to calculate land surface temperature"""

    # test different input types
    surf_temp = Meteorology.calc_surf_temp(
        dew_point=params["dew_point"], h_swe=params["h_swe"]
    )

    assert np.allclose(surf_temp, -2.5)

    # test different input data
    surf_temp = Meteorology.calc_surf_temp(dew_point=[-2, 2], h_swe=[0, 1])

    assert_almost_equal(surf_temp, [-2, 0])


@pytest.mark.parametrize(
    "params",
    [
        PARAMS["test_case_1"],
        PARAMS["test_case_2"],
        PARAMS["test_case_3"],
        PARAMS["test_case_4"],
    ],
)
def test_calc_precipitable_water_content(params):
    """Test to calculate precipitable water content"""

    prec_water_content = Meteorology.calc_precipitable_water_content(
        dew_point=params["dew_point"]
    )

    assert np.allclose(prec_water_content, 0.96062486)


@pytest.mark.parametrize(
    "params",
    [
        PARAMS["test_case_1"],
        PARAMS["test_case_2"],
        PARAMS["test_case_3"],
        PARAMS["test_case_4"],
    ],
)
def test_calc_latent_heat_flux(params):
    """Test to calculate latent heat flux"""

    latent_heat_flux = Meteorology.calc_latent_heat_flux(
        rho_air=1.2614,
        aero_conductance=params["aero_conductance"],
        air_pressure=params["air_pressure"],
        air_vapor_pressure=params["air_vapor_pressure"],
        surf_vapor_pressure=params["surf_vapor_pressure"],
    )

    assert np.allclose(latent_heat_flux, 9.09061757)


def test_calc_tsn_offset():
    """Test to calculate the true solar noon offset"""

    tsn_offset = Meteorology.calc_tsn_offset(
        julian_day=2.5, year=2023, longitude=[-75, 75], gmt_offset=[-5, 5]
    )
    assert np.allclose(tsn_offset, 0.07271961)


# TODO: check shortwave radiation
@pytest.mark.parametrize(
    "params",
    [
        PARAMS["test_case_1"],
        # PARAMS['test_case_2'],
        # PARAMS['test_case_3'],
        # PARAMS['test_case_4']
    ],
)
def test_calc_net_shortwave_radiation(params):
    """Test to calculate net shortwave energy flux"""

    # clear sky
    q_sw = Meteorology.calc_net_shortwave_radiation(
        julian_day=2.5,
        tsn_offset=params["tsn_offset"],
        latitude=params["latitude"],
        prec_water=params["prec_water"],
        alpha=params["alpha"],
        beta=params["beta"],
        albedo=params["albedo"],
        canopy_factor=params["canopy_factor"],
        cloud_factor=params["cloud_factor"],
        dust_atten=params["dust_atten"],
    )

    assert np.allclose(q_sw, 257.9950632421464)

    # consider cloud and canopy impact
    q_sw = Meteorology.calc_net_shortwave_radiation(
        julian_day=2.5,
        tsn_offset=params["tsn_offset"],
        latitude=params["latitude"],
        prec_water=params["prec_water"],
        alpha=params["alpha"],
        beta=params["beta"],
        albedo=params["albedo"],
        dust_atten=params["dust_atten"],
        canopy_factor=params["canopy_factor"],
        cloud_factor=params["cloud_factor"],
        clear_sky=False,
    )

    assert np.allclose(q_sw, 3.42945865680005)


@pytest.mark.parametrize(
    "params",
    [
        PARAMS["test_case_1"],
        PARAMS["test_case_2"],
        PARAMS["test_case_3"],
        PARAMS["test_case_4"],
    ],
)
def test_calc_air_emissivity(params):
    """Test to calculate emissivity of air"""

    # Brutsaert method
    air_em = Meteorology.calc_air_emissivity(
        air_temp=params["air_temp"],
        air_vapor_pressure=params["air_vapor_pressure"],
        canopy_factor=params["canopy_factor"],
        cloud_factor=params["cloud_factor"],
        satterlund=False,
    )
    assert np.allclose(air_em, 0.99526025)

    # Satterlund method
    air_em = Meteorology.calc_air_emissivity(
        air_temp=params["air_temp"],
        air_vapor_pressure=params["air_vapor_pressure"],
        canopy_factor=params["canopy_factor"],
        cloud_factor=params["cloud_factor"],
        satterlund=True,
    )
    assert np.allclose(air_em, 0.77507235)


@pytest.mark.parametrize(
    "params",
    [
        PARAMS["test_case_1"],
        PARAMS["test_case_2"],
        PARAMS["test_case_3"],
        PARAMS["test_case_4"],
    ],
)
def test_calc_net_longwave_radiation(params):
    """Test to calculate net long wave energy flux"""

    q_lw = Meteorology.calc_net_longwave_radiation(
        air_temp=params["air_temp"],
        surf_temp=params["surf_temp"],
        air_emissivity=params["air_emissivity"],
        surf_emissivity=params["surf_emissivity"],
    )
    assert np.allclose(q_lw, 17.69477686)


@pytest.mark.parametrize(
    "q_lw, q_sw, qh, qe",
    [
        (1, 2, 3, 4),
        ([1], [2], [3], [4]),
        ([1, 1], [2, 2], [3, 3], [4, 4]),
        ([[1, 1], [1, 1]], [[2, 2], [2, 2]], [[3, 3], [3, 3]], [[4, 4], [4, 4]]),
    ],
)
def test_calc_net_energy_flux(q_lw, q_sw, qh, qe):
    """Test to calculate total net energy flux"""

    q_sum = Meteorology.calc_net_energy_flux(
        shortwave_energy_flux=q_sw,
        longwave_energy_flux=q_lw,
        sensible_heat_flux=qh,
        latent_heat_flux=qe,
    )

    assert np.allclose(q_sum, 10)


def test_run_multiple_steps():
    """Test to run multiple steps"""

    # step 1: run one step with default setting
    grid = RasterModelGrid((2, 2))
    grid.add_full("atmosphere_bottom_air__temperature", 8.0, at="node")
    grid.add_full("land_surface__temperature", 4.0, at="node")
    grid.add_full("land_surface__latitude", 40, at="node")
    grid.add_full("land_surface__longitude", -75, at="node")
    grid.add_full("atmosphere_bottom_air__pressure", 1000, at="node")
    grid.add_full("atmosphere_bottom_air_flow__reference-height_speed", 2.0)
    grid.add_full("atmosphere_bottom_air_water-vapor__relative_saturation", 0.5)
    grid.add_full("atmosphere_bottom_air__brutsaert_emissivity_canopy_factor", 0.98)
    grid.add_full("atmosphere_bottom_air__brutsaert_emissivity_cloud_factor", 0.62)
    grid.add_full("atmosphere_aerosol_dust__reduction_of_transmittance", 0.1)
    grid.add_full("land_surface__slope_angle", 0.2)
    grid.add_full("land_surface__aspect_angle", 0.5)
    grid.add_full("land_surface__albedo", 0.3)

    dt = 60 * 60  # 1hr

    met = Meteorology(grid, start_datetime="2023-01-02 12:00:00", gmt_offset=-5)
    met.run_one_step(dt)

    assert_almost_equal(
        grid.at_node["land_surface_net-shortwave-radiation__energy_flux"], 3.34685251
    )
    assert_almost_equal(
        grid.at_node["land_surface_net-longwave-radiation__energy_flux"], 17.6928257
    )
    assert_almost_equal(
        grid.at_node["atmosphere_bottom_air_land_net-sensible-heat__energy_flux"],
        9.743864,
    )
    assert_almost_equal(
        grid.at_node["atmosphere_bottom_air_land_net-latent-heat__energy_flux"],
        4.9008869,
    )
    assert_almost_equal(
        grid.at_node["land_surface_net-total-energy__energy_flux"], 35.6844291
    )  # Q_sum

    # step2: change input values
    grid.at_node["atmosphere_bottom_air__temperature"].fill(-1)
    grid.at_node["land_surface__temperature"].fill(2)
    grid.at_node["land_surface__albedo"].fill(0.1)
    grid.at_node["atmosphere_bottom_air__brutsaert_emissivity_canopy_factor"].fill(0.1)
    grid.at_node["atmosphere_bottom_air__brutsaert_emissivity_cloud_factor"].fill(0.2)

    met.run_one_step(dt)

    assert_almost_equal(
        grid.at_node["land_surface_net-shortwave-radiation__energy_flux"], 166.10733404
    )
    assert_almost_equal(
        grid.at_node["land_surface_net-longwave-radiation__energy_flux"], -109.48594478
    )
    assert_almost_equal(
        grid.at_node["atmosphere_bottom_air_land_net-sensible-heat__energy_flux"],
        -121.49949326,
    )
    assert_almost_equal(
        grid.at_node["atmosphere_bottom_air_land_net-latent-heat__energy_flux"],
        -43.26555593,
    )
    assert_almost_equal(
        grid.at_node["land_surface_net-total-energy__energy_flux"], -108.14365992
    )

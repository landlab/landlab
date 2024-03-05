import numpy as np

from numpy.testing import (
    assert_almost_equal,
    assert_array_almost_equal,
    assert_array_equal,
    assert_equal,
    assert_raises,
)

from landlab import HexModelGrid, RadialModelGrid, RasterModelGrid
from landlab.components.genveg.photosynthesis import Photosynthesis


def test_calculate_hourly_direct_light_extinction():
    """Test direct light attenuation coefficient calculation"""
    solar_elevation = 0.467396
    kdr_teh = 1.1097223
    kdr_genveg = Photosynthesis.calculate_hourly_direct_light_extinction(
        solar_elevation
    )
    assert_almost_equal(kdr_genveg, kdr_teh)


def test_calculate_hourly_diffuse_light_extinction():
    """Test diffuse light attenuation coeffcient calculation"""
    lai = 0.0244633
    kdf_teh = 0.96219749
    kdf_genveg = Photosynthesis.calculate_hourly_diffuse_light_extinction(lai)
    assert_almost_equal(kdf_genveg, kdf_teh)

def test_calculate_incremental_PAR():
    """Test Gaussian integration of PAR into hourly increments"""
    increment_hour = 
    solar_elevation = 
    grid_par_W_per_sqm =  
    current_day = 
    incremental_direct_PAR_teh = 
    incremental_diffuse_PAR_teh = 
    hourly_direct_PAR, hourly_diffuse_PAR = Photosynthesis.calculate_incremental_PAR(
    increment_hour, solar_elevation, grid_par_W_per_sqm, current_day
    )
    assert_almost_equal(hourly_direct_PAR, incremental_direct_PAR_teh)
    assert_almost_equal(hourly_diffuse_PAR, incremental_diffuse_PAR_teh)


def test_absorbed_incremental_par():
    #    """Test part of photosynthesis algorithm is producing output similar to Teh method
    #    AbsorbedHourPAR in solar.h lines 64-87
    #
    #
    #    """
    #    increment_hour = 17.7097
    #
    #    solar_elevation = 0.467396
    #    grid_par_W_per_sqm =
    lai = 0.012


#    _current_day =
#    absorbed_hour_PAR_sunlit = 80.0121489
#    absorbed_hour_PAR_shaded = 24.7014599

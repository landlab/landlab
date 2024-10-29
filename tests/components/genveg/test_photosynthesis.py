import numpy as np
from numpy.testing import assert_allclose

from landlab.components.genveg.photosynthesis import Photosynthesis

photo_object = Photosynthesis(0.9074, _current_day=90)

dtypes = [
    ("species", "U10"),
    ("pid", int),
    ("cell_index", int),
    ("x_loc", float),
    ("y_loc", float),
    (("root", "root_biomass"), float),
    (("leaf", "leaf_biomass"), float),
    (("stem", "stem_biomass"), float),
    (("reproductive", "repro_biomass"), float),
    ("dead_root", float),
    ("dead_stem", float),
    ("dead_leaf", float),
    ("dead_reproductive", float),
    ("dead_age", float),
    ("shoot_sys_width", float),
    ("root_sys_width", float),
    ("shoot_sys_height", float),
    ("root_sys_depth", float),
    ("total_leaf_area", float),
    ("live_leaf_area", float),
    ("plant_age", float),
    ("n_stems", int),
    ("pup_x_loc", float),
    ("pup_y_loc", float),
    ("pup_cost", float),
    ("item_id", int),
]
plants = np.empty((0, 26), dtype=dtypes)
plantlist = []
plant = "Corn"
pidval = 0
cell_index = 0
for i in range(1):
    pidval = i
    cell_index = i + 1
    plantlist.append(
        (
            plant,
            pidval,
            cell_index,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0,
            np.nan,
            np.nan,
            np.nan,
            i,
        )
    )
plants = np.array(plantlist, dtype=dtypes)
plants["root"] = np.array([0.80000000000000004])
plants["stem"] = np.array([0.29999999999999999])
plants["leaf"] = np.array([0.50000000000000000])
plants["live_leaf_area"] = 0.022 * plants["leaf"]
plants["shoot_sys_width"] = np.array(
    [(4 * (plants["live_leaf_area"] / 0.012) / np.pi) ** 0.5]
)
plants["shoot_sys_height"] = np.array([0.010000000000000000])
plants["total_leaf_area"] = plants["live_leaf_area"]
plants["plant_age"] = np.array([90])
plants["n_stems"] = np.array([1.0])


def test_calculate_hourly_direct_light_extinction():
    """Test direct light attenuation coefficient calculation"""
    solar_elevation = np.array([0.467396])
    kdr_teh = np.array([1.1097223])
    kdr_genveg = photo_object.calculate_hourly_direct_light_extinction(solar_elevation)
    assert_allclose(kdr_genveg, kdr_teh, rtol=0.0001)


def test_calculate_hourly_diffuse_light_extinction():
    """Test diffuse light attenuation coeffcient calculation"""
    lai = np.array([0.0244633, 0.012])
    kdf_teh = np.array([0.96219749, 0.97307909472236043])
    kdf_genveg = photo_object.calculate_hourly_diffuse_light_extinction(lai)
    assert_allclose(kdf_genveg, kdf_teh, rtol=0.0001)


def test_calculate_hour_temp():
    increment_hour = np.array([12.00, 18.461916000000002])
    min_temp = np.array([0.10000000000000001])
    max_temp = np.array([16.699999999999999])
    hour_temp_teh = np.array([15.552812015345607, 6.0881570941648402])
    hour_temp_calc_noon = photo_object.calculate_hour_temp(
        increment_hour[0], min_temp, max_temp
    )
    hour_temp_calc_afternoon = photo_object.calculate_hour_temp(
        increment_hour[1], min_temp, max_temp
    )
    assert_allclose(hour_temp_calc_noon, hour_temp_teh[0], rtol=0.0001)
    assert_allclose(hour_temp_calc_afternoon, hour_temp_teh[1], rtol=0.0001)


def test_calculate_incremental_PAR():
    """Test Gaussian integration of PAR into hourly increments"""
    increment_hour = np.array([12.00, 18.461916000000002])
    solar_elevation = np.array([0.72481400133962604, -0.025788369553778639])
    grid_par_W_per_sqm = np.array(
        [18399966.731363025 / 86400, 18399966.731363025 / 86400]
    )
    current_day = 90  # march 31
    incremental_direct_PAR_teh = np.array([444.04291756 / 2, 0.00])
    incremental_diffuse_PAR_teh = np.array([196.4907107 / 2, 0.00])
    hourly_direct_PAR, hourly_diffuse_PAR = photo_object.calculate_incremental_PAR(
        increment_hour, solar_elevation, grid_par_W_per_sqm, current_day
    )
    assert_allclose(hourly_direct_PAR, incremental_direct_PAR_teh, rtol=0.0001)
    assert_allclose(hourly_diffuse_PAR, incremental_diffuse_PAR_teh, rtol=0.0001)


def test_absorbed_incremental_par():
    """
    Test part of photosynthesis algorithm is producing output similar to Teh method
    AbsorbedHourPAR in solar.h lines 64-87
    """
    increment_hour = np.array([6.2902924124321347, 8.6070323903780466])
    solar_elevation = np.array([0.095170832314239964, 0.45119603341512748])
    grid_par_W_per_sqm = np.array(
        [18399966.731363025 / 86400 / 2, 18399966.731363025 / 86400 / 2]
    )
    lai = 0.012
    current_day = 90
    absorbed_hour_PAR_sunlit_teh = np.array(
        [82.012148983212555 * 4.55, 182.23112770728699 * 4.55]
    )
    absorbed_hour_PAR_shaded_teh = np.array(
        [24.701459944977820 * 4.55, 55.227013122938175 * 4.55]
    )
    absorbed_hour_PAR_sunlit, absorbed_hour_PAR_shaded = (
        photo_object.calculate_absorbed_incremental_PAR(
            increment_hour, solar_elevation, grid_par_W_per_sqm, lai, current_day
        )
    )
    assert_allclose(absorbed_hour_PAR_shaded, absorbed_hour_PAR_shaded_teh, rtol=0.0001)
    assert_allclose(absorbed_hour_PAR_sunlit, absorbed_hour_PAR_sunlit_teh, rtol=0.0001)


def test_calculate_leaf_assimilation():
    """
    Test leaf assimilation calculation in assim.h from Teh with assumption that leaf
    temperature is the same as air temperature
    """
    CO2_conc = np.array([245])
    hour_temp = photo_object.calculate_hour_temp(np.array([12]))
    sunlit_par = np.array([951.36361799661461])
    shaded_par = np.array([341.88946400925738])
    shaded_leaf_assim_teh = np.array([10.100904909580670])
    sunlit_leaf_assim_teh = np.array([28.107427842696641])
    sunlit_assim = photo_object.calculate_leaf_assimilation(
        sunlit_par, CO2_conc, hour_temp
    )
    shaded_assim = photo_object.calculate_leaf_assimilation(
        shaded_par, CO2_conc, hour_temp
    )
    assert_allclose(sunlit_assim, sunlit_leaf_assim_teh, rtol=0.0001)
    assert_allclose(shaded_assim, shaded_leaf_assim_teh, rtol=0.0001)


def test_calculate_sunlit_shaded_lai_proportion():
    """
    Test the subdivision of leaf area index to sunlit and shaded based on
    solar elevation
    """
    solar_elevation = np.array(
        [0.095170832314239964, 0.45119603341512748, 0.72481400133962604]
    )
    lai = np.array([0.012, 0.012, 0.012])
    sunlit_lai_teh = np.array(
        [0.011629010198938363, 0.011917816564890738, 0.011945864475940510]
    )
    shaded_lai_teh = np.array(
        [0.00037098980106163755, 8.2183435109262418e-05, 5.4135524059490542e-05]
    )
    sunlit_lai, shaded_lai = photo_object.calculate_sunlit_shaded_lai_proportion(
        solar_elevation, lai
    )
    assert_allclose(sunlit_lai, sunlit_lai_teh, rtol=0.0001)
    assert_allclose(shaded_lai, shaded_lai_teh, rtol=0.0001)


def test_photosynthesize():
    """
    Test the integrated output of GenVeg photosythesize with the equivalent output
    from Teh modified to output per plant

    NEED TO AMEND FOR STOMATAL CONDUCTANCE ADJUSTMENT
    """
    grid_par_W_per_sqm = np.array([18399966.731363025 / 86400 / 2])
    min_temp = np.array([0.10000000000000001])
    max_temp = np.array([16.699999999999999])

    lai = np.array([0.012])
    gphot = photo_object.photosynthesize(
        grid_par_W_per_sqm, min_temp, max_temp, lai, plants, 90
    )
    gphot_per_unit_area_teh = np.array([0.24187992803952649 / 0.60266080364221775])
    gphot_plant_teh = gphot_per_unit_area_teh * (
        np.pi / 4 * plants["shoot_sys_width"] ** 2
    )
    assert_allclose(gphot, gphot_plant_teh, rtol=0.01)

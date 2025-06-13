import numpy as np
import pytest
from numpy.testing import assert_allclose
from numpy.testing import assert_array_almost_equal
from numpy.testing import assert_equal

from landlab.components.genveg.allometry import Biomass
from landlab.components.genveg.habit import Habit
from landlab.components.genveg.species import Species

dt = np.timedelta64(1, 'D')


def create_species_object(example_input_params):
    return Species(example_input_params["BTS"], dt=dt, latitude=0.9074)


def create_habit_object(example_input_params):
    params = example_input_params["BTS"]
    params["grow_params"]["abg_biomass"] = {}
    params["morph_params"]["canopy_area"] = {}
    params["grow_params"]["abg_biomass"]["min"] = params["grow_params"]["plant_part_min"]["leaf"] + params["grow_params"]["plant_part_min"]["stem"]
    params["grow_params"]["abg_biomass"]["max"] = params["grow_params"]["plant_part_max"]["leaf"] + params["grow_params"]["plant_part_max"]["stem"]
    params["morph_params"]["canopy_area"]["min"] = params["morph_params"]["shoot_sys_width"]["min"]**2 * 0.25 * np.pi
    params["morph_params"]["canopy_area"]["max"] = params["morph_params"]["shoot_sys_width"]["max"]**2 * 0.25 * np.pi
    allometry = Biomass(params, empirical_coeffs={"root_dia_coeffs": {"a": 0.08, "b": 0.24}})
    return Habit(params, allometry, dt=1, green_parts=("leaf", "stem"))


def test_calc_canopy_area_from_shoot_width(example_input_params):
    h = create_habit_object(example_input_params)
    # zero array returns zero array
    assert_array_almost_equal(
        h._calc_canopy_area_from_shoot_width(shoot_sys_width=np.array([0, 0, 0])),
        np.array([0, 0, 0]),
    )
    # single input returns correct input
    assert_equal(
        h._calc_canopy_area_from_shoot_width(shoot_sys_width=0.325), 0.08295768100885548
    )
    # an array of values
    assert_allclose(
        h._calc_canopy_area_from_shoot_width(
            shoot_sys_width=np.array([0, 0.0004, 0.678, 1.5, 3])
        ),
        np.array(
            [
                0.00000000e00,
                1.25663706e-07,
                3.61034969e-01,
                1.76714587e00,
                7.06858347e00,
            ]
        ),
    )


def test__calc_diameter_from_area(example_input_params):
    h = create_habit_object(example_input_params)
    canopy_area = np.array([0.1, 1.28, 3.7])
    shoot_width = np.array([0.35682, 1.27662, 2.17048])
    pred_shoot_width = h._calc_diameter_from_area(canopy_area)
    assert_array_almost_equal(shoot_width, pred_shoot_width, decimal=5)


def test__calc_canopy_volume(example_input_params):
    h = create_habit_object(example_input_params)
    shoot_width = np.array([0.08, 0.35, 0.72])
    height = np.array([0.1, 0.5, 1.7])
    basal_dia = np.array([0.005, 0.02, 0.1])
    volume = np.array([0.00017868, 0.017004, 0.2672])
    pred_vol = h._calc_canopy_volume(shoot_width, basal_dia, height)
    assert_allclose(pred_vol, volume, rtol=0.0001)


def test_calc_root_sys_width(example_input_params):
    h = create_habit_object(example_input_params)
    shoot_width = np.array([0.08, 0.35, 0.72])
    height = np.array([0.1, 0.5, 1.7])
    basal_dia = np.array([0.005, 0.02, 0.1])
    root_width = np.array([0.080043, 0.084081, 0.144128])
    pred_root_width = h.calc_root_sys_width(shoot_width, basal_dia, height)
    assert_allclose(pred_root_width, root_width, rtol=0.0001)


def test_estimate_abg_biomass_from_cover(example_input_params, example_plant_array):
    h = create_habit_object(example_input_params)
    abg_biomass = example_plant_array["leaf"] + example_plant_array["stem"]
    example_plant_array["basal_dia"], example_plant_array["shoot_sys_width"], example_plant_array["shoot_sys_height"] = h.calc_abg_dims_from_biomass(abg_biomass)
    est_abg_biomass = h.estimate_abg_biomass_from_cover(example_plant_array)
    assert_allclose(est_abg_biomass, abg_biomass, rtol=0.0001)


def test_calc_canopy_area_from_shoot_width_raises_error(example_input_params):
    h = create_habit_object(example_input_params)
    with pytest.raises(ValueError):
        h._calc_canopy_area_from_shoot_width(-1.5)
        h._calc_canopy_area_from_shoot_width(np.array([0, 0.004, -0.678, 1.5, 3]))

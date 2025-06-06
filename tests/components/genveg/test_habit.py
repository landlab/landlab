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
    allometry = Biomass(params)
    return Habit(params, allometry, dt=1, green_parts=("leaf", "stem"))

"""
def test_calc_canopy_area_from_shoot_width(example_input_params):
    h = Habit(params=example_input_params["BTS"], green_parts=("leaf", "stem"))
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
"""


def test__calc_diameter_from_area(example_input_params):
    h = create_habit_object(example_input_params)
    canopy_area = np.array([0.1, 1.28, 3.7])
    shoot_width = np.array([0.35682, 1.27662, 2.17048])
    pred_shoot_width = h._calc_diameter_from_area(canopy_area)
    assert_array_almost_equal(shoot_width, pred_shoot_width, decimal=5)

"""
def test__calc_abg_dims_from_biomass(example_input_params):
    h = Habit(params=example_input_params["BTS"], green_parts=("leaf", "stem"))
    abg_biomass = np.array([2.1, 0.2, 15.2])
    


def test_calc_crown_area_from_shoot_width_raises_error(example_input_params):
    create_species_object(example_input_params)
    h = Habit(params=example_input_params["BTS"], green_parts=("leaf", "stem"))
    with pytest.raises(ValueError):
        h._calc_canopy_area_from_shoot_width(-1.5)
        h._calc_canopy_area_from_shoot_width(np.array([0, 0.004, -0.678, 1.5, 3]))
"""
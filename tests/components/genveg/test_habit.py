import numpy as np
import pytest
from numpy.testing import assert_allclose
from numpy.testing import assert_array_almost_equal
from numpy.testing import assert_equal

from landlab.components.genveg.habit import Habit
from landlab.components.genveg.species import Species

dt = np.timedelta64(1, 'D')


def create_species_object(example_input_params):
    return Species(example_input_params["BTS"], dt=dt, latitude=0.9074)


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


def test__calc2_allometry_coeffs(example_input_params):
    h = Habit(params=example_input_params["BTS"], green_parts=("leaf"))
    x_min = np.array([1.])
    x_max = np.array([3.])
    y_min = np.array([2.])
    y_max = np.array([6.])
    slope = np.array([1])
    y_intercept = np.array([np.log(2)])
    (b, m) = h._calc2_allometry_coeffs(x_min, x_max, y_min, y_max)
    assert_array_almost_equal(y_intercept, b, decimal=5)
    assert_array_almost_equal(slope, m, decimal=5)


def test__calc3_allometry_coeffs(example_input_params):
    h = Habit(params=example_input_params["BTS"], green_parts=("leaf"))
    x_min = np.array([1.1])
    x_mean = np.array([1.5])
    x_max = np.array([3.])
    y_min = np.array([2.])
    y_mean = np.array([4.])
    y_max = np.array([6.])
    z_min = np.array([5.2])
    z_mean = np.array([9.5])
    z_max = np.array([16.7])
    # Assuming form is ln(z) = a + b ln(x) + c ln(y)
    a_wa = np.array([1.13487])
    b_wa = np.array([0.413507])
    c_wa = np.array([0.684388])
    (a, b, c) = h._calc3_allometry_coeffs(x_min, x_mean, x_max, y_min, y_mean, y_max, z_min, z_mean, z_max)
    assert_array_almost_equal(a, a_wa, decimal=5)
    assert_array_almost_equal(b, b_wa, decimal=5)
    assert_array_almost_equal(c, c_wa, decimal=5)


def test_apply2_allometry_eq_for_xs(example_input_params):
    h = Habit(params=example_input_params["BTS"], green_parts=("leaf"))
    h.morph_params["basal_coeffs"] = {"a": np.log(2), "b": 1}
    abg_biomass = np.array([2.1, 0.2, 15.2])
    pred_basal_diameters = h._apply2_allometry_eq_for_xs(abg_biomass, "basal_coeffs")
    basal_diameters = np.array([1.05, 0.1, 7.6])
    assert_array_almost_equal(pred_basal_diameters, basal_diameters, decimal=5)


def test_apply2_allometry_eq_for_ys(example_input_params):
    h = Habit(params=example_input_params["BTS"], green_parts=("leaf"))
    h.morph_params["basal_coeffs"] = {"a": np.log(2), "b": 1}
    basal_diameter = np.array([0.5, 3.5, 12])
    pred_abg_biomass = h._apply2_allometry_eq_for_ys(basal_diameter, "basal_coeffs")
    abg_biomass = np.array([1., 7., 24.])
    assert_array_almost_equal(pred_abg_biomass, abg_biomass, decimal=5)


def test_apply3_allometry_eq_for_ys(example_input_params):
    h = Habit(params=example_input_params["BTS"], green_parts=("leaf", "stem"))
    h.morph_params["canopy_coeffs"] = {"a": 1.13487, "b": 0.413507, "c": 0.684388}
    abg_biomass = np.array([2.1, 0.2, 15.2])
    basal_dia = np.array([0.5, 3.5, 12])
    pred_canopy = h._apply3_allometry_eq_for_ys(basal_dia, abg_biomass, "canopy_coeffs")
    canopy = np.array([0.856126, 0.008508, 2.262885])
    assert_array_almost_equal(pred_canopy, canopy, decimal=6)


def test_apply3_allometry_eq_for_xs(example_input_params):
    h = Habit(params=example_input_params["BTS"], green_parts=("leaf", "stem"))
    h.morph_params["canopy_coeffs"] = {"a": 1.13487, "b": 0.413507, "c": 0.684388}
    abg_biomass = np.array([2.1, 0.2, 15.2])
    canopy_area = np.array([1.28, 0.05, 3.7])
    pred_basal_dia = h._apply3_allometry_eq_for_xs(canopy_area, abg_biomass, "canopy_coeffs")
    basal_dia = np.array([0.256964, 0.186657, 5.318099])
    assert_array_almost_equal(pred_basal_dia, basal_dia, decimal=6)


def test_apply3_allometry_eq_for_zs(example_input_params):
    h = Habit(params=example_input_params["BTS"], green_parts=("leaf", "stem"))
    h.morph_params["canopy_coeffs"] = {"a": 1.13487, "b": 0.413507, "c": 0.684388}
    canopy_area = np.array([1.28, 0.05, 3.7])
    basal_dia = np.array([0.5, 3.5, 12])
    pred_abg = h._apply3_allometry_eq_for_zs(basal_dia, canopy_area, "canopy_coeffs")
    abg = np.array([2.765432, 0.672101, 21.280764])
    assert_array_almost_equal(pred_abg, abg, decimal=6)


def test__calc_diameter_from_area(example_input_params):
    h = Habit(params=example_input_params["BTS"], green_parts=("leaf", "stem"))
    canopy_area = np.array([0.1, 1.28, 3.7])
    shoot_width = np.array([0.35682, 1.27662, 2.17048])
    pred_shoot_width = h._calc_diameter_from_area(canopy_area)
    assert_array_almost_equal(shoot_width, pred_shoot_width, decimal=5)


def test__calc_abg_dims_from_biomass(example_input_params):
    h = Habit(params=example_input_params["BTS"], green_parts=("leaf", "stem"))
    abg_biomass = np.array([2.1, 0.2, 15.2])
    


def test_calc_crown_area_from_shoot_width_raises_error(example_input_params):
    create_species_object(example_input_params)
    h = Habit(params=example_input_params["BTS"], green_parts=("leaf", "stem"))
    with pytest.raises(ValueError):
        h._calc_canopy_area_from_shoot_width(-1.5)
        h._calc_canopy_area_from_shoot_width(np.array([0, 0.004, -0.678, 1.5, 3]))

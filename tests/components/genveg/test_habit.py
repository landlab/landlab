import numpy as np
import pytest
from numpy.testing import assert_array_almost_equal, assert_equal, assert_allclose

from landlab.components.genveg.habit import Habit
from landlab.components.genveg.species import Species


def test_calc_canopy_area_from_shoot_width(example_input_params):
    Species(
        species_params=example_input_params["BTS"],
        latitude=0.9074
    )
    h = Habit(
        params=example_input_params['BTS'],
        green_parts=('leaf', 'stem')
    )
    # zero array returns zero array
    assert_array_almost_equal(
        h._calc_canopy_area_from_shoot_width(shoot_sys_width=np.array([0, 0, 0])),
        np.array([0, 0, 0])
    )
    # single input returns correct input
    assert_equal(
        h._calc_canopy_area_from_shoot_width(shoot_sys_width=0.325),
        0.08295768100885548
    )
    # an array of values
    assert_allclose(
        h._calc_canopy_area_from_shoot_width(shoot_sys_width=np.array([0, .0004, 0.678, 1.5, 3])),
        np.array([0.00000000e+00, 1.25663706e-07, 3.61034969e-01, 1.76714587e+00, 7.06858347e+00])
    )

import numpy as np
import pytest

from numpy.testing import assert_allclose
from landlab.components.genveg.species import Species


def create_species_object(example_input_params):
    return Species(species_params=example_input_params["BTS"], latitude=0.9074)


def test_get_daily_nsc_concentration(example_input_params):
    species_object = create_species_object(example_input_params)
    days = example_input_params["BTS"]["duration_params"]
    parts = ["root", "leaf", "reproductive", "stem"]
    nsc_content_gs_start_actual = {
        "root": 276.9606782 / 1000,
        "leaf": 247.3937921 / 1000,
        "stem": 100 / 1000,
        "reproductive": 287.459626 / 1000,
    }
    nsc_content_spring195_actual = {
        "root": 205.575914 / 1000,
        "leaf": 223.822083 / 1000,
        "stem": 93.9498113 / 1000,
        "reproductive": 220.297893 / 1000,
    }
    gs_start_day = days["growing_season_start"]
    nsc_content_gs_start = species_object.get_daily_nsc_concentration(gs_start_day)
    nsc_content_195 = species_object.get_daily_nsc_concentration(195)
    for part in parts:
        assert_allclose(
            nsc_content_gs_start_actual[part], nsc_content_gs_start[part], rtol=0.0001
        )
        assert_allclose(
            nsc_content_spring195_actual[part], nsc_content_195[part], rtol=0.0001
        )


def test_caulate_dervied_params_error_message_raise(example_input_params):
    species_object = create_species_object(example_input_params)
    example_input_params["BTS"]["morph_params"]["max_shoot_sys_width"] = np.negative(
        example_input_params["BTS"]["morph_params"]["max_shoot_sys_width"]
    )
    with pytest.raises(ValueError):
        species_object.calculate_derived_params(example_input_params["BTS"])


import numpy as np

from numpy.testing import assert_allclose, assert_almost_equal

from landlab.components.genveg.species import Species
from landlab.components.genveg.habit import Forbherb, Graminoid, Shrub, Tree, Vine
from landlab.components.genveg.form import Bunch, Colonizing, Multiplestems, Rhizomatous, Singlecrown, Singlestem, Stoloniferous, Thicketforming


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


def test_calc_area_of_circle(example_input_params):
    species_object = create_species_object(example_input_params)
    morph_params = example_input_params["BTS"]["morph_params"]
    m_params = ["max_shoot_sys_width", "min_shoot_sys_width", "max_root_sys_width", "min_root_sys_width"]
    # values from excel sheet
    area_values = np.array([0.070685835, 0.0000785398, 0.096211275, 0.0000785398])
    for m_param, a_value in zip(m_params, area_values):
        assert_almost_equal(
            species_object.calc_area_of_circle(morph_params[m_param]),
            a_value
        )


# Test calculate_derived_params functions
def test_max_vitial_volume(example_input_params):
    assert_almost_equal(
        create_species_object(example_input_params).calc_volume_cylinder(
            area=0.070685835,
            height=example_input_params["BTS"]["morph_params"]["max_height"]
        ),
        0.053014376
    )


def test_ratio_calculations(example_input_params):
    species_object = create_species_object(example_input_params)
    morph_param = example_input_params["BTS"]["morph_params"]
    # area_per_stem
    assert_almost_equal(
        species_object.calc_param_ratio(0.070685835, morph_param["max_n_stems"]),
        0.007068583
    )
    # min_abg_aspect_ratio
    assert_almost_equal(
        species_object.calc_param_ratio(morph_param["max_height"], morph_param["min_shoot_sys_width"]),
        75
    )
    # max_abg_aspect_ratio
    assert_almost_equal(
        species_object.calc_param_ratio(morph_param["max_height"], morph_param["max_shoot_sys_width"]),
        2.5
    )
    # min_basal_ratio
    assert_almost_equal(
        species_object.calc_param_ratio(morph_param["min_shoot_sys_width"], morph_param["min_basal_dia"]),
        1.843906736
    )
    # max_basal_ratio
    assert_almost_equal(
        species_object.calc_param_ratio(morph_param["max_shoot_sys_width"], morph_param["max_basal_dia"]),
        3
    )
    # biomass_packing
    assert_almost_equal(
        species_object.calc_param_ratio(17.9, 0.053014376),
        337.6442646
    )
    # senesce_rate
    assert_almost_equal(
        species_object.calc_param_ratio(0.9, 32),
        0.028125
    )


def test_sum_vars_in_calculate_derived_params(example_input_params):
    species_object = create_species_object(example_input_params)
    species_param = species_object.calculate_derived_params(example_input_params["BTS"])
    # Checked via excel
    # Max total Biomass
    assert_almost_equal(
        species_param["grow_params"]["max_total_biomass"],
        17.9
    )
    # max_growth_biomass
    assert_almost_equal(
        species_param["grow_params"]["max_growth_biomass"],
        13.9
    )
    # max_abg_biomass
    assert_almost_equal(
        species_param["grow_params"]["max_abg_biomass"],
        9.6
    )
    # min_total_biomass
    assert_almost_equal(
        species_param["grow_params"]["min_total_biomass"],
        0.062222222
    )
    # min_growth_biomass
    assert_almost_equal(
        species_param["grow_params"]["min_growth_biomass"],
        0.062222222
    )
    # min_abg_biomass
    assert_almost_equal(
        species_param["grow_params"]["min_abg_biomass"],
        0.052222222
    )
    # min_nsc_biomass
    assert_almost_equal(
        species_param["grow_params"]["min_nsc_biomass"],
        0.03369
    )


def test_nsc_rate_change_per_season_and_part(example_input_params):
    species_object = create_species_object(example_input_params)
    species_param = species_object.calculate_derived_params(example_input_params["BTS"])
    ncs_rate_change = species_param["duration_params"]["nsc_rate_change"]
    # winter_nsc_rate
    # - leaf
    assert_almost_equal(
        ncs_rate_change["winter_nsc_rate"]["leaf"],
        0.003676471
    )
    # - reproductive
    assert_almost_equal(
        ncs_rate_change["winter_nsc_rate"]["reproductive"],
        -0.004595588
    )
    # - root
    assert_almost_equal(
        ncs_rate_change["winter_nsc_rate"]["root"],
        -0.003676471
    )
    # - stem
    assert_almost_equal(
        ncs_rate_change["winter_nsc_rate"]["stem"],
        -0.00245098
    )
    # spring_nsc_rate
    # - leaf
    assert_almost_equal(
        ncs_rate_change["spring_nsc_rate"]["leaf"],
        -0.015060241
    )
    # - reproductive
    assert_almost_equal(
        ncs_rate_change["spring_nsc_rate"]["reproductive"],
        -0.041415663
    )
    # - root
    assert_almost_equal(
        ncs_rate_change["spring_nsc_rate"]["root"],
        -0.045180723
    )
    # - stem
    assert_almost_equal(
        ncs_rate_change["spring_nsc_rate"]["stem"],
        -0.006024096
    )
    # summer_nsc_rate
    # - leaf
    assert_almost_equal(
        ncs_rate_change["summer_nsc_rate"]["leaf"],
        -0.02173913
    )
    # - reproductive
    assert_almost_equal(
        ncs_rate_change["summer_nsc_rate"]["reproductive"],
        0.042119565
    )
    # - root
    assert_almost_equal(
        ncs_rate_change["summer_nsc_rate"]["root"],
        0.054347826
    )
    # - stem
    assert_almost_equal(
        ncs_rate_change["summer_nsc_rate"]["stem"],
        0.010869565
    )
    # fall_nsc_rate
    # - leaf
    assert_almost_equal(
        ncs_rate_change["fall_nsc_rate"]["leaf"],
        0.046875
    )
    # - reproductive
    assert_almost_equal(
        ncs_rate_change["fall_nsc_rate"]["reproductive"],
        0.076171875
    )
    # - root
    assert_almost_equal(
        ncs_rate_change["fall_nsc_rate"]["root"],
        0.0625
    )
    # - stem
    assert_almost_equal(
        ncs_rate_change["fall_nsc_rate"]["stem"],
        0.015625
    )


def test_select_habit_class(example_input_params):
    species_object = create_species_object(example_input_params)
    dummy_species = example_input_params
    for spec, cls in zip(
        ['forb_herb', 'graminoid', 'shrub', 'tree', 'vine'],
        [Forbherb, Graminoid, Shrub, Tree, Vine]
    ):
        dummy_species["BTS"]["plant_factors"]["growth_habit"] = spec
        print(dummy_species["BTS"]["plant_factors"]["growth_habit"])
        assert isinstance(
            species_object.select_habit_class(dummy_species["BTS"]),
            cls
        )



def test_select_form_class(example_input_params):
    species_object = create_species_object(example_input_params)
    dummy_growth_form = example_input_params
    for growth, cls in zip(
        ['bunch', 'colonizing', 'multiple_stems', 'rhizomatous', 'single_crown', 'single_stem', 'stoloniferous', 'thicket_forming'],
        [Bunch, Colonizing, Multiplestems, Rhizomatous, Singlecrown, Singlestem, Stoloniferous, Thicketforming]
    ):
        dummy_growth_form["BTS"]["plant_factors"]["growth_form"] = growth
        assert isinstance(
            species_object.select_form_class(dummy_growth_form["BTS"]),
            cls
        )

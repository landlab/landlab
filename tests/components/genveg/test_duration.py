import numpy as np
from numpy.testing import assert_almost_equal, assert_array_less

from landlab.components.genveg.duration import Duration, Annual, Perennial, Evergreen, Deciduous


dt = np.timedelta64(1, 'D')


def create_duration_object(example_input_params, green_parts=("leaf", "stem") , persistent_parts=("root"), child=None):
    grow_params = example_input_params["BTS"]["grow_params"]
    duration_params = example_input_params["BTS"]["duration_params"]
    child_class = {
        "annual": Annual(grow_params, duration_params, dt=dt),
        "perennial": Perennial(grow_params, duration_params, dt, green_parts, persistent_parts),
        "evergreen": Evergreen(grow_params, duration_params, dt),
        "deciduous": Deciduous(grow_params, duration_params, dt, green_parts),
    }
    if child is not None:
        return child_class[child]
    else:
        return Duration(grow_params, duration_params, dt=dt, green_parts=green_parts, persistent_parts=persistent_parts)


def test_enter_dormancy(example_input_params, example_plant):
    dur_object = create_duration_object(example_input_params)
    initial_leaf = example_plant["leaf"].copy()
    initial_root = example_plant["root"].copy()
    plant_out = dur_object.enter_dormancy(example_plant)
    assert_almost_equal(plant_out["leaf"], np.zeros_like(initial_leaf), decimal=6)
    assert_almost_equal(plant_out["root"], initial_root, decimal=6)


def test_emerge(example_input_params, example_plant):
    # Annual
    annual_object = create_duration_object(example_input_params, child="annual")
    example_plant["root"] = example_plant["leaf"] = example_plant["stem"] = np.array([0.])
    max_start = 2 * annual_object.growdict["growth_min_biomass"]
    min_start = annual_object.growdict["growth_min_biomass"]
    available_stored_biomass = total_persistent_biomass = np.zeros_like(example_plant["root"])
    emerged_plant = annual_object.emerge(example_plant, available_stored_biomass, total_persistent_biomass)
    emerged_plant_total = emerged_plant["root"] + emerged_plant["stem"] + emerged_plant["leaf"]
    assert_array_less(min_start, emerged_plant_total)
    assert_array_less(emerged_plant_total, max_start)
    # Evergreen
    evergreen_object = create_duration_object(example_input_params, child="evergreen")
    emerged_plant = evergreen_object.emerge(example_plant, available_stored_biomass, total_persistent_biomass)
    assert_almost_equal(emerged_plant["root"], example_plant["root"])
    assert_almost_equal(emerged_plant["leaf"], example_plant["leaf"])
    assert_almost_equal(emerged_plant["stem"], example_plant["stem"])
    # Deciduous
    decid_object = create_duration_object(example_input_params, child="deciduous")
    example_plant["root"] = 0.8
    available_stored_biomass = example_plant["root"] * 0.33333
    total_persistent_biomass = example_plant["root"]
    emerged_plant = decid_object.emerge(example_plant, available_stored_biomass, total_persistent_biomass)
    mass_green = emerged_plant["stem"] + emerged_plant["leaf"]
    assert_array_less(mass_green, available_stored_biomass)


def test_set_initial_biomass(example_input_params, example_plant):
    # Testing for Deciduous set_initial_biomass, testing for other classes handled under test_emerge
    decid_object = create_duration_object(example_input_params, child="deciduous")
    in_growing_season = True
    example_plant["root"] = example_plant["leaf"] = example_plant["stem"] = np.array([0.])
    max_start = 2 * decid_object.growdict["growth_min_biomass"]
    min_start = decid_object.growdict["growth_min_biomass"]
    emerged_plant = decid_object.set_initial_biomass(example_plant, in_growing_season)
    emerged_plant_total = emerged_plant["root"] + emerged_plant["stem"] + emerged_plant["leaf"]
    assert_array_less(min_start, emerged_plant_total)
    assert_array_less(emerged_plant_total, max_start)
    in_growing_season = False
    min_root = decid_object.growdict["plant_part_min"]["root"]
    max_root = 2 * min_root
    emerged_plant = decid_object.set_initial_biomass(example_plant, in_growing_season)
    assert_array_less(min_root, emerged_plant["root"])
    assert_array_less(emerged_plant["root"], max_root)
    assert_almost_equal(emerged_plant["leaf"], np.zeros_like(emerged_plant["leaf"]))
    assert_almost_equal(emerged_plant["stem"], np.zeros_like(emerged_plant["stem"]))


def test__solve_biomass_allocation(example_input_params, example_plant):
    dur_object = create_duration_object(example_input_params)
    total_biomass = example_plant["root"] + example_plant["leaf"] + example_plant["stem"]
    # Values from Excel spreadsheet applying Poorter formulations
    ideal_root = np.array([0.562235])
    ideal_leaf = np.array([0.621115])
    ideal_stem = np.array([0.416640])
    root, leaf, stem = dur_object._solve_biomass_allocation(total_biomass)
    assert_almost_equal(root, ideal_root, decimal=5)
    assert_almost_equal(leaf, ideal_leaf, decimal=5)
    assert_almost_equal(stem, ideal_stem, decimal=5)

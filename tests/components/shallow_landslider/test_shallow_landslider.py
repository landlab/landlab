# tests/components/test_shallow_landslider.py

import numpy as np

from landlab import RasterModelGrid
from landlab.components import ShallowLandslider


def make_grid(ny=5, nx=5, spacing=10.0):
    mg = RasterModelGrid((ny, nx), xy_spacing=spacing)
    z = mg.add_zeros("topographic__elevation", at="node")
    z[:] = np.linspace(0, ny * nx - 1, ny * nx)
    mg.add_full("soil__depth", 0.5, at="node")
    mg.add_full("earthquake__horizontal_pga", 0.3, at="node")
    mg.add_full("earthquake__vertical_pga", 0.1, at="node")
    return mg


def test_initialization_only_creates_output_field():
    mg = make_grid()
    fields_before = set(mg.at_node)

    ShallowLandslider(mg, cohesion_eff=100, angle_int_frict=30)

    assert set(mg.at_node) == fields_before | {"landslide__selected_labels"}


def test_run_one_step_pipeline():
    mg = make_grid()

    comp = ShallowLandslider(
        mg,
        cohesion_eff=25.0,
        angle_int_frict=27.0,
    )

    comp.run_one_step(
        kde_input={
            "kde_data": {"overall": None},
            "kde_transform": {"x_bounds": (0, 1)},
            "width_threshold": 1.5,
        }
    )

    sel = mg.at_node["landslide__selected_labels"]
    asp = comp.results["aspect_labels"]
    assert set(np.unique(sel)) <= set(np.unique(asp))


def test_run_one_step_without_measured_width_data():
    mg = make_grid()

    comp = ShallowLandslider(
        mg,
        cohesion_eff=25.0,
        angle_int_frict=27.0,
    )

    comp.run_one_step()

    assert comp.results["split_labels"] is None
    assert "landslide__dimension_split_labels" not in mg.at_node
    assert "landslide__selected_labels" in mg.at_node


def test_probabilistic_selection_sets_a_proportion():
    mg = make_grid()
    comp = ShallowLandslider(mg, cohesion_eff=25.0, angle_int_frict=27.0)
    comp.run_one_step()
    assert comp._selected_proportion is not None


def test_results_property_contains_expected_keys():
    mg = make_grid()
    comp = ShallowLandslider(mg, cohesion_eff=25.0, angle_int_frict=27.0)
    comp.run_one_step()

    r = comp.results
    for key in [
        "factor_of_safety",
        "a_transient",
        "a_driving",
        "a_diff",
        "unstable_mask",
        "labels",
        "aspect_labels",
        "selected_labels",
        "group_properties",
    ]:
        assert key in r


def test_aspect_labels_refine_region_labels():
    mg = make_grid()
    comp = ShallowLandslider(mg, cohesion_eff=15, angle_int_frict=30)
    comp.run_one_step()

    reg = comp.results["labels"]
    asp = comp.results["aspect_labels"]

    # Aspect subgrouping should not introduce labels with no parent region.
    assert set(np.unique(asp)) - {0} <= set(np.unique(reg)) - {0}


def test_dimension_split_labels_consistency():
    mg = make_grid()
    comp = ShallowLandslider(
        mg,
        cohesion_eff=20,
        angle_int_frict=30,
    )
    comp.run_one_step(
        kde_input={
            "kde_data": {"overall": None},
            "kde_transform": {},
            "width_threshold": 2.0,
        }
    )

    asp = comp.results["aspect_labels"]
    dim = comp.results["split_labels"]

    # Split labels should refine but never contradict aspect groups
    assert set(np.unique(dim)) - {0} <= set(np.unique(asp)) - {0}

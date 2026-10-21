# tests/helper_functions/test_regions_selection_split.py

import numpy as np
import pytest

from landlab import RasterModelGrid
from landlab.components import ShallowLandslider
from landlab.components.shallow_landslider.shallow_landslide_component import (
    _regionprops,
)


def make_grid(shape=(6, 6), spacing=10.0):
    mg = RasterModelGrid(shape, xy_spacing=spacing)
    z = mg.add_zeros("topographic__elevation", at="node")
    z[:] = np.arange(z.size)
    soil = mg.add_zeros("soil__depth", at="node")
    soil[:] = 1.0
    mg.add_zeros("earthquake__horizontal_pga", at="node")
    mg.add_zeros("earthquake__vertical_pga", at="node")
    return mg


def test_regionprops_matches_expected_geometry():
    labels = np.zeros((5, 6), dtype=int)
    labels[1:3, 2:5] = 2

    (region,) = _regionprops(labels)

    assert region.label == 2
    assert region.bbox == (1, 2, 3, 5)
    assert region.area == 6.0
    assert np.allclose(region.centroid, (1.5, 3.0))
    assert np.isclose(region.axis_major_length, 4.0 * np.sqrt(2.0 / 3.0))
    assert np.isclose(region.axis_minor_length, 2.0)
    assert np.isclose(region.orientation, np.pi / 2.0)
    assert np.isclose(region.eccentricity, np.sqrt(5.0 / 8.0))


def test_calculate_regions_connectivity(monkeypatch):
    mg = make_grid()
    comp = ShallowLandslider(mg, cohesion_eff=10, angle_int_frict=30)

    mask = np.zeros(mg.shape, dtype=bool)
    mask[1, 1] = True
    mask[2, 2] = True

    labels4, n4 = comp._calculate_regions(mask, connect_val=4)
    labels8, n8 = comp._calculate_regions(mask, connect_val=8)

    assert n4 == 2
    assert n8 == 1


def test_fill_region_holes_preserves_raw_labels():
    mg = make_grid(shape=(7, 7))
    comp = ShallowLandslider(mg, cohesion_eff=10, angle_int_frict=30)
    labels = np.zeros(mg.shape, dtype=int)
    labels[1:6, 1:6] = 1
    labels[2:5, 2:5] = 0
    comp._labels = labels.ravel()

    comp._fill_region_holes()

    assert np.array_equal(comp._labels.reshape(mg.shape), labels)
    assert np.all(comp._filled_labels.reshape(mg.shape)[1:6, 1:6] == 1)
    assert np.count_nonzero(comp._hole_fill_mask) == 9


def test_fill_region_holes_retains_open_cavities():
    mg = make_grid(shape=(7, 7))
    comp = ShallowLandslider(mg, cohesion_eff=10, angle_int_frict=30)
    labels = np.zeros(mg.shape, dtype=int)
    labels[1:6, 1:6] = 1
    labels[2:5, 2:5] = 0
    labels[1, 3] = 0
    comp._labels = labels.ravel()

    comp._fill_region_holes()

    assert np.array_equal(comp._filled_labels.reshape(mg.shape), labels)
    assert not np.any(comp._hole_fill_mask)


@pytest.mark.parametrize("excluded", ["nodata", "closed_node", "other_region"])
def test_fill_region_holes_preserves_excluded_cells(excluded):
    mg = make_grid(shape=(7, 7))
    comp = ShallowLandslider(mg, cohesion_eff=10, angle_int_frict=30)
    labels = np.zeros(mg.shape, dtype=int)
    labels[1:6, 1:6] = 1
    labels[2:5, 2:5] = 0
    if excluded == "nodata":
        nodata = mg.add_zeros("nodata__mask", at="node", dtype=bool)
        nodata.reshape(mg.shape)[3, 3] = True
    elif excluded == "closed_node":
        mg.status_at_node[3 * mg.shape[1] + 3] = mg.BC_NODE_IS_CLOSED
    else:
        labels[3, 3] = 2
    comp._labels = labels.ravel()

    comp._fill_region_holes()

    filled = comp._filled_labels.reshape(mg.shape)
    assert filled[3, 3] == labels[3, 3]
    surrounding = np.ones((3, 3), dtype=bool)
    surrounding[1, 1] = False
    assert np.all(filled[2:5, 2:5][surrounding] == 1)
    assert np.count_nonzero(comp._hole_fill_mask) == 8


def test_zone_split_by_aspect():
    mg = make_grid()
    comp = ShallowLandslider(mg, cohesion_eff=10, angle_int_frict=30)

    groups = np.zeros(mg.shape, dtype=int)
    groups[1:5, 1:5] = 1

    aspect = np.zeros(mg.shape)
    aspect[:3, :] = 10.0
    aspect[3:, :] = 190.0

    zones = comp._create_zones(90)
    subgroups, zone_labels, info = comp._split_groups_by_aspect(
        groups, aspect, zones=zones, min_size=2
    )
    assert subgroups.max() >= 2
    assert np.any(subgroups[:3, :] > 0) and np.any(subgroups[3:, :] > 0)


def test_probabilistic_group_selection_reproducible():
    mg = make_grid()
    comp = ShallowLandslider(mg, cohesion_eff=10, angle_int_frict=30)

    labeled = np.zeros(mg.shape, dtype=int)
    labeled[1:3, 1:3] = 1
    labeled[3:5, 3:5] = 2

    probs = np.zeros_like(labeled, dtype=float)
    probs[labeled == 1] = 0.8
    probs[labeled == 2] = 0.2

    sel1, meta1 = comp._probabilistic_group_selection(
        labeled, probs, random_seed=123, reproducible=True
    )
    sel2, meta2 = comp._probabilistic_group_selection(
        labeled, probs, random_seed=123, reproducible=True
    )

    assert np.array_equal(sel1, sel2)
    assert meta1["proportion_calculated"] == meta2["proportion_calculated"]


class FakeKDE:
    def resample(self, n):
        return np.vstack([np.ones(n) * 5.0, np.ones(n) * 10.0])


def test_recursive_split_wide_regions_splits():
    mg = make_grid()
    comp = ShallowLandslider(mg, cohesion_eff=10, angle_int_frict=30)

    labeled = np.zeros(mg.shape, dtype=int)
    labeled[1:5, 1:5] = 1

    aspect = np.ones(mg.shape) * 45
    slopes = np.ones(mg.shape) * 20
    kde_results = {"overall": FakeKDE()}
    transform_info = {"log_x": False, "log_y": False}

    new_labels, info = comp._recursive_split_wide_regions(
        labeled,
        aspect,
        slopes,
        kde_results,
        transform_info,
        width_threshold=0.5,
        max_iterations=1,
        min_region_size=5,
        convergence_threshold=0.95,
    )

    assert new_labels.max() > 1


def test_recursive_split_converges_early():
    mg = make_grid()
    comp = ShallowLandslider(mg, cohesion_eff=10, angle_int_frict=30)

    labels = np.zeros(mg.shape, dtype=int)
    labels[1:5, 1:5] = 1

    aspect = np.ones(mg.shape) * 45
    slopes = np.ones(mg.shape) * 20

    class TinyKDE:
        def resample(self, n):
            # produce widths matching actual so no split
            return np.vstack([np.ones(n) * 100.0, np.ones(n) * 100.0])

    kde = {"overall": TinyKDE()}
    info = {"log_x": False, "log_y": False}

    new_labels, splits = comp._recursive_split_wide_regions(
        labels,
        aspect,
        slopes,
        kde,
        info,
        width_threshold=2.0,
        max_iterations=5,
        min_region_size=3,
        convergence_threshold=0.9,
    )
    assert splits == []  # no splits → early convergence


def test_probabilistic_group_selection_all_zero_probabilities():
    mg = make_grid()
    comp = ShallowLandslider(mg, cohesion_eff=10, angle_int_frict=30)

    labels = np.zeros(mg.shape, dtype=int)
    labels[1:3, 1:3] = 1
    probs = np.zeros_like(labels, dtype=float)

    sel, meta = comp._probabilistic_group_selection(labels, probs, reproducible=True)
    assert meta["proportion_calculated"] == 0.0
    assert np.all(sel == 0)

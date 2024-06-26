#!/usr/bin/env python
"""Tests for SpeciesEvolver zone objects."""
import numpy as np
import pandas as pd
import pytest

from landlab import RasterModelGrid
from landlab.components import SpeciesEvolver
from landlab.components.species_evolution import ZoneController
from landlab.components.species_evolution import ZoneTaxon
from landlab.components.species_evolution import zone as zn


@pytest.fixture()
def zone_example_grid():
    mg = RasterModelGrid((5, 7), 2)
    z = mg.add_zeros("node", "topographic__elevation")
    return mg, z


def zone_func(grid):
    z = grid.at_node["topographic__elevation"]
    return z == 1


def zone_func_with_vars(grid, var1, var2):
    z = grid.at_node["topographic__elevation"]
    return np.all([z == 1, grid.x_of_node > var1, grid.y_of_node > var2], 0)


def test_none_to_none(zone_example_grid):
    mg, z = zone_example_grid

    sc = ZoneController(mg, zone_func)
    sc.run_one_step(1)

    np.testing.assert_array_equal(len(sc.zones), 0)


def test_none_to_one(zone_example_grid):
    mg, z = zone_example_grid

    # No zones exist at time 0.

    se = SpeciesEvolver(mg)
    sc = ZoneController(mg, zone_func)
    taxa = sc.populate_zones_uniformly(1)
    se.track_taxa(taxa)

    np.testing.assert_equal(len(sc.zones), 0)

    expected_df = pd.DataFrame(
        {
            "time": [0],
            "zones": [0],
            "fragmentations": [np.nan],
            "captures": [np.nan],
            "area_captured_sum": [np.nan],
            "area_captured_max": [np.nan],
        }
    )
    pd.testing.assert_frame_equal(sc.record_data_frame, expected_df, check_like=True)

    # Create a zone for time 1.

    z[[9, 10, 11, 12]] = 1

    sc.run_one_step(1)
    se.run_one_step(1)

    np.testing.assert_equal(len(sc.zones), 1)
    np.testing.assert_equal(sc.zones[0]._conn_type, zn.Connection.NONE_TO_ONE)

    expected_df = pd.DataFrame(
        {
            "time": [0, 1],
            "zones": [0, 1],
            "fragmentations": [np.nan, 0],
            "captures": [np.nan, 0],
            "area_captured_sum": [np.nan, 0],
            "area_captured_max": [np.nan, 0],
        }
    )
    pd.testing.assert_frame_equal(sc.record_data_frame, expected_df, check_like=True)

    np.testing.assert_equal(se.record_data_frame.taxa.sum(), 0)


def test_one_to_none(zone_example_grid):
    mg, z = zone_example_grid

    # Create a zone for time 0.

    z[[9, 10, 11, 12]] = 1

    se = SpeciesEvolver(mg)
    sc = ZoneController(mg, zone_func)
    taxa = sc.populate_zones_uniformly(1)
    se.track_taxa(taxa)

    np.testing.assert_equal(len(sc.zones), 1)
    zone = sc.zones[0]

    expected_df = pd.DataFrame(
        {
            "time": [0],
            "zones": [1],
            "fragmentations": [np.nan],
            "captures": [np.nan],
            "area_captured_sum": [np.nan],
            "area_captured_max": [np.nan],
        }
    )
    pd.testing.assert_frame_equal(sc.record_data_frame, expected_df, check_like=True)

    # No zones for time 1.

    z[:] = 0

    sc.run_one_step(1)
    se.run_one_step(1)

    np.testing.assert_equal(len(sc.zones), 0)
    np.testing.assert_equal(zone._conn_type, zn.Connection.ONE_TO_NONE)

    expected_df = pd.DataFrame(
        {
            "time": [0, 1],
            "zones": [1, 0],
            "fragmentations": [np.nan, 0],
            "captures": [np.nan, 0],
            "area_captured_sum": [np.nan, 0],
            "area_captured_max": [np.nan, 0],
        }
    )
    pd.testing.assert_frame_equal(sc.record_data_frame, expected_df, check_like=True)

    np.testing.assert_equal(se.record_data_frame.taxa.sum(), 1)


def test_one_to_one(zone_example_grid):
    mg, z = zone_example_grid

    # Create a zone for time 0.

    z[[9, 10, 11, 12]] = 1

    se = SpeciesEvolver(mg)
    sc = ZoneController(mg, zone_func)
    taxa = sc.populate_zones_uniformly(1)
    se.track_taxa(taxa)

    np.testing.assert_equal(len(sc.zones), 1)
    zone = sc.zones[0]

    expected_df = pd.DataFrame(
        {
            "time": [0],
            "zones": [1],
            "fragmentations": [np.nan],
            "captures": [np.nan],
            "area_captured_sum": [np.nan],
            "area_captured_max": [np.nan],
        }
    )
    pd.testing.assert_frame_equal(sc.record_data_frame, expected_df, check_like=True)

    # Modify elevation, although  there is still one zone in time 1.

    z[[11, 12]] = 0

    sc.run_one_step(1)
    se.run_one_step(1)

    np.testing.assert_equal(len(sc.zones), 1)
    np.testing.assert_equal(zone, sc.zones[0])
    np.testing.assert_equal(sc.zones[0]._conn_type, zn.Connection.ONE_TO_ONE)

    expected_df = pd.DataFrame(
        {
            "time": [0, 1],
            "zones": [1, 1],
            "fragmentations": [np.nan, 0],
            "captures": [np.nan, 0],
            "area_captured_sum": [np.nan, 0],
            "area_captured_max": [np.nan, 0],
        }
    )
    pd.testing.assert_frame_equal(sc.record_data_frame, expected_df, check_like=True)

    np.testing.assert_equal(len(se.get_extant_taxon_objects(time=1)), 1)


def test_one_to_many(zone_example_grid):
    mg, z = zone_example_grid

    # Create a zone for time 0.

    z[[9, 10, 11, 12]] = 1

    se = SpeciesEvolver(mg)
    sc = ZoneController(mg, zone_func)
    taxa = sc.populate_zones_uniformly(1)
    se.track_taxa(taxa)

    expected_df = pd.DataFrame(
        {
            "time": [0],
            "zones": [1],
            "fragmentations": [np.nan],
            "captures": [np.nan],
            "area_captured_sum": [np.nan],
            "area_captured_max": [np.nan],
        }
    )
    pd.testing.assert_frame_equal(sc.record_data_frame, expected_df, check_like=True)

    np.testing.assert_equal(len(se.get_extant_taxon_objects(time=0)), 1)

    # Break the zone in two for time 1.

    z[11] = 0

    sc.run_one_step(1)
    se.run_one_step(1)

    np.testing.assert_equal(len(sc.zones), 2)
    np.testing.assert_equal(
        {z._conn_type for z in sc.zones}, {None, zn.Connection.ONE_TO_MANY}
    )

    expected_df = pd.DataFrame(
        {
            "time": [0, 1],
            "zones": [1, 2],
            "fragmentations": [np.nan, 2],
            "captures": [np.nan, 0],
            "area_captured_sum": [np.nan, 0],
            "area_captured_max": [np.nan, 0],
        }
    )
    pd.testing.assert_frame_equal(sc.record_data_frame, expected_df, check_like=True)

    np.testing.assert_equal(len(se.get_extant_taxon_objects()), 2)


def test_many_to_one(zone_example_grid):
    mg, z = zone_example_grid

    # Create two zones for time 0.

    z[[8, 9, 10, 12]] = 1

    se = SpeciesEvolver(mg)
    sc = ZoneController(mg, zone_func)
    taxa = sc.populate_zones_uniformly(1)
    se.track_taxa(taxa)

    expected_df = pd.DataFrame(
        {
            "time": [0],
            "zones": [2],
            "fragmentations": [np.nan],
            "captures": [np.nan],
            "area_captured_sum": [np.nan],
            "area_captured_max": [np.nan],
        }
    )
    pd.testing.assert_frame_equal(sc.record_data_frame, expected_df, check_like=True)

    np.testing.assert_equal(len(se.get_extant_taxon_objects(time=0)), 2)

    # Modify elevation such that two zones each overlap the original two zones.

    z[11] = 1

    sc.run_one_step(1)
    se.run_one_step(1)

    np.testing.assert_equal(len(sc.zones), 1)
    np.testing.assert_equal(sc.zones[0]._conn_type, zn.Connection.MANY_TO_ONE)

    expected_df = pd.DataFrame(
        {
            "time": [0, 1],
            "zones": [2, 1],
            "fragmentations": [np.nan, 0],
            "captures": [np.nan, 1],
            "area_captured_sum": [np.nan, 12],
            "area_captured_max": [np.nan, 12],
        }
    )
    pd.testing.assert_frame_equal(sc.record_data_frame, expected_df, check_like=True)

    np.testing.assert_equal(len(se.get_extant_taxon_objects(time=1)), 2)


def test_many_to_many(zone_example_grid):
    mg, z = zone_example_grid

    # Create two zones for time 0.

    z[[10, 12, 17, 19, 24, 26]] = 1

    se = SpeciesEvolver(mg)
    sc = ZoneController(mg, zone_func)
    taxa = sc.populate_zones_uniformly(1)
    se.track_taxa(taxa)

    expected_df = pd.DataFrame(
        {
            "time": [0],
            "zones": [2],
            "fragmentations": [np.nan],
            "captures": [np.nan],
            "area_captured_sum": [np.nan],
            "area_captured_max": [np.nan],
        }
    )
    pd.testing.assert_frame_equal(sc.record_data_frame, expected_df, check_like=True)

    np.testing.assert_equal(len(se.get_extant_taxon_objects(time=0)), 2)

    # Modify elevation such that two zones each overlap the original two zones.

    z[[17, 19]] = 0
    z[[11, 25]] = 1

    sc.run_one_step(1)
    se.run_one_step(1)

    np.testing.assert_equal(len(sc.zones), 2)
    for zone in sc.zones:
        np.testing.assert_equal(zone._conn_type, zn.Connection.MANY_TO_MANY)

    expected_df = pd.DataFrame(
        {
            "time": [0, 1],
            "zones": [2, 2],
            "fragmentations": [np.nan, 0],
            "captures": [np.nan, 2],
            "area_captured_sum": [np.nan, 24],
            "area_captured_max": [np.nan, 12],
        }
    )
    pd.testing.assert_frame_equal(sc.record_data_frame, expected_df, check_like=True)

    np.testing.assert_equal(len(se.get_extant_taxon_objects()), 4)


def test_one_to_many_to_one(zone_example_grid):
    mg, z = zone_example_grid

    z[[9, 10, 11, 12]] = 1

    se = SpeciesEvolver(mg)
    sc = ZoneController(mg, zone_func)
    taxa = sc.populate_zones_uniformly(1, time_to_allopatric_speciation=1)
    se.track_taxa(taxa)

    z[11] = 0
    sc.run_one_step(1)
    se.run_one_step(1)

    z[11] = 1
    sc.run_one_step(1)
    se.run_one_step(1)

    np.testing.assert_equal(len(se.get_extant_taxon_objects()), 1)


def test_min_area(zone_example_grid):
    mg, z = zone_example_grid

    # Create a zone for time 0.

    z[[9, 11, 12]] = 1

    sc = ZoneController(mg, zone_func, minimum_area=5)
    np.testing.assert_equal(len(sc.zones), 1)


def test_neighborhood_structure(zone_example_grid):
    mg, z = zone_example_grid

    z[[10, 16]] = 1

    sc = ZoneController(mg, zone_func)
    np.testing.assert_equal(len(sc.zones), 1)

    sc = ZoneController(mg, zone_func, neighborhood_structure="D4")
    np.testing.assert_equal(len(sc.zones), 2)

    np.testing.assert_raises(
        ValueError, ZoneController, mg, zone_func, neighborhood_structure="D"
    )


def test_zone_func_kwargs(zone_example_grid):
    mg, z = zone_example_grid
    z[mg.core_nodes] = 1
    sc = ZoneController(mg, zone_func_with_vars, var1=1, var2=2)

    expected_mask = np.array(
        [
            False,
            False,
            False,
            False,
            False,
            False,
            False,
            False,
            False,
            False,
            False,
            False,
            False,
            False,
            False,
            True,
            True,
            True,
            True,
            True,
            False,
            False,
            True,
            True,
            True,
            True,
            True,
            False,
            False,
            False,
            False,
            False,
            False,
            False,
            False,
        ]
    )
    np.testing.assert_array_equal(sc.zones[0].mask, expected_mask)

    z[-2:] = 0
    expected_mask[-2:] = False

    sc.run_one_step(10)

    np.testing.assert_array_equal(sc.zones[0].mask, expected_mask)


def test_zone_taxon_range_mask(zone_example_grid):
    mg, z = zone_example_grid

    # Create a zone for time 0.

    ids = [9, 10, 11, 12]
    z[ids] = 1

    se = SpeciesEvolver(mg)
    sc = ZoneController(mg, zone_func)
    taxa = sc.populate_zones_uniformly(1)
    se.track_taxa(taxa)

    expected_mask = np.zeros(mg.number_of_nodes, bool)
    expected_mask[ids] = True
    np.testing.assert_array_equal(taxa[0].range_mask, expected_mask)

    # Remove extent so taxa range mask is all False.

    z[ids] = 0
    sc.run_one_step(1)
    se.run_one_step(1)

    expected_mask = np.zeros(mg.number_of_nodes, bool)
    np.testing.assert_array_equal(taxa[0].range_mask, expected_mask)


def test_time_to_allopatric_speciation(zone_example_grid):
    mg, z = zone_example_grid

    # Create a zone for time 0.

    z[[9, 10, 11, 12]] = 1

    se = SpeciesEvolver(mg)
    sc = ZoneController(mg, zone_func)
    taxa = sc.populate_zones_uniformly(1, time_to_allopatric_speciation=20)
    se.track_taxa(taxa)

    z[[11]] = 0

    while len(se.taxa_data_frame) == 1:
        sc.run_one_step(10)
        se.run_one_step(10)

    expected_df = pd.DataFrame(
        {
            "pid": [np.nan, 0],
            "type": 2 * [ZoneTaxon.__name__],
            "t_first": [0, 30],
            "t_final": 2 * [np.nan],
        },
        index=[0, 1],
    )
    expected_df.index.name = "tid"
    expected_df["pid"] = expected_df["pid"].astype("Int64")
    expected_df["t_final"] = expected_df["t_final"].astype("Int64")

    pd.testing.assert_frame_equal(se.taxa_data_frame, expected_df, check_like=True)


def test_pseudoextinction(zone_example_grid):
    mg, z = zone_example_grid

    z[[9, 10, 11, 12]] = 1

    se = SpeciesEvolver(mg)
    sc = ZoneController(mg, zone_func)
    taxa = sc.populate_zones_uniformly(1, persists_post_speciation=False)
    se.track_taxa(taxa)

    z[11] = 0
    sc.run_one_step(1)
    se.run_one_step(1)

    expected_df = pd.DataFrame(
        {
            "pid": [np.nan, 0, 0],
            "type": 3 * [ZoneTaxon.__name__],
            "t_first": [0, 1, 1],
            "t_final": [1, np.nan, np.nan],
        },
        index=[0, 1, 2],
    )
    expected_df.index.name = "tid"
    expected_df["pid"] = expected_df["pid"].astype("Int64")
    expected_df["t_final"] = expected_df["t_final"].astype("Int64")

    pd.testing.assert_frame_equal(se.taxa_data_frame, expected_df, check_like=True)

    expected_df = pd.DataFrame(
        {
            "time": [0, 1],
            "taxa": [1, 2],
            "speciations": [np.nan, 2],
            "extinctions": [np.nan, 0],
            "pseudoextinctions": [np.nan, 1],
        }
    )
    pd.testing.assert_frame_equal(se.record_data_frame, expected_df, check_like=True)

#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""Tests for SpeciesEvolver zone objects."""
import numpy as np
import pandas as pd
import pytest

from landlab import RasterModelGrid
from landlab.components import SpeciesEvolver
from landlab.components.species_evolution import ZoneController
from landlab.components.species_evolution import zone as zn


@pytest.fixture()
def zone_example_grid():
    mg = RasterModelGrid((3, 7), 2)
    z = mg.add_zeros('node', 'topographic__elevation')
    return mg, z


def zone_func(grid):
    z = grid.at_node['topographic__elevation']
    return z == 1


def zone_func_with_vars(grid, var1, var2):
    z = grid.at_node['topographic__elevation']
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
    species = sc.populate_zones_uniformly(1)
    se.introduce_species(species)

    np.testing.assert_equal(len(sc.zones), 0)

    expected_df = pd.DataFrame({
        'time': [0],
        'zone_count': [0],
        'fragmentation_count': [np.nan],
        'capture_count': [np.nan],
        'area_captured_sum': [np.nan],
        'area_captured_max': [np.nan]}
    )
    pd.testing.assert_frame_equal(
        sc.record_data_frame, expected_df, check_like=True
    )

    # Create a zone for time 1.

    z[[9, 10, 11, 12]] = 1

    sc.run_one_step(1)
    se.run_one_step(1)

    np.testing.assert_equal(len(sc.zones), 1)
    np.testing.assert_equal(sc.zones[0]._conn_type, zn._NONE_TO_ONE)

    expected_df = pd.DataFrame({
        'time': [0, 1],
        'zone_count': [0, 1],
        'fragmentation_count': [np.nan, 0],
        'capture_count': [np.nan, 0],
        'area_captured_sum': [np.nan, 0],
        'area_captured_max': [np.nan, 0]}
    )
    pd.testing.assert_frame_equal(
        sc.record_data_frame, expected_df, check_like=True
    )

    np.testing.assert_equal(se.record_data_frame.species_count.sum(), 0)


def test_one_to_none(zone_example_grid):
    mg, z = zone_example_grid

    # Create a zone for time 0.

    z[[9, 10, 11, 12]] = 1

    se = SpeciesEvolver(mg)
    sc = ZoneController(mg, zone_func)
    species = sc.populate_zones_uniformly(1)
    se.introduce_species(species)

    np.testing.assert_equal(len(sc.zones), 1)
    zone = sc.zones[0]

    expected_df = pd.DataFrame({
        'time': [0],
        'zone_count': [1],
        'fragmentation_count': [np.nan],
        'capture_count': [np.nan],
        'area_captured_sum': [np.nan],
        'area_captured_max': [np.nan]}
    )
    pd.testing.assert_frame_equal(
        sc.record_data_frame, expected_df, check_like=True
    )

    # No zones for time 1.

    z[:] = 0

    sc.run_one_step(1)
    se.run_one_step(1)

    np.testing.assert_equal(len(sc.zones), 0)
    np.testing.assert_equal(zone._conn_type, zn._ONE_TO_NONE)

    expected_df = pd.DataFrame({
        'time': [0, 1],
        'zone_count': [1, 0],
        'fragmentation_count': [np.nan, 0],
        'capture_count': [np.nan, 0],
        'area_captured_sum': [np.nan, 0],
        'area_captured_max': [np.nan, 0]}
    )
    pd.testing.assert_frame_equal(
        sc.record_data_frame, expected_df, check_like=True
    )

    np.testing.assert_equal(se.record_data_frame.species_count.sum(), 1)


def test_one_to_one(zone_example_grid):
    mg, z = zone_example_grid

    # Create a zone for time 0.

    z[[9, 10, 11, 12]] = 1

    se = SpeciesEvolver(mg)
    sc = ZoneController(mg, zone_func)
    species = sc.populate_zones_uniformly(1)
    se.introduce_species(species)

    np.testing.assert_equal(len(sc.zones), 1)
    zone = sc.zones[0]

    expected_df = pd.DataFrame({
        'time': [0],
        'zone_count': [1],
        'fragmentation_count': [np.nan],
        'capture_count': [np.nan],
        'area_captured_sum': [np.nan],
        'area_captured_max': [np.nan]}
    )
    pd.testing.assert_frame_equal(
        sc.record_data_frame, expected_df, check_like=True
    )

    # Modify elevation, although  there is still one zone in time 1.

    z[[11, 12]] = 0

    sc.run_one_step(1)
    se.run_one_step(1)

    np.testing.assert_equal(len(sc.zones), 1)
    np.testing.assert_equal(zone, sc.zones[0])
    np.testing.assert_equal(sc.zones[0]._conn_type, zn._ONE_TO_ONE)

    expected_df = pd.DataFrame({
        'time': [0, 1],
        'zone_count': [1, 1],
        'fragmentation_count': [np.nan, 0],
        'capture_count': [np.nan, 0],
        'area_captured_sum': [np.nan, 0],
        'area_captured_max': [np.nan, 0]}
    )
    pd.testing.assert_frame_equal(
        sc.record_data_frame, expected_df, check_like=True
    )

    np.testing.assert_equal(len(se.filter_species(time=1)), 1)


def test_one_to_many(zone_example_grid):
    mg, z = zone_example_grid

    # Create a zone for time 0.

    z[[9, 10, 11, 12]] = 1

    se = SpeciesEvolver(mg)
    sc = ZoneController(mg, zone_func)
    species = sc.populate_zones_uniformly(1)
    se.introduce_species(species)

    expected_df = pd.DataFrame({
        'time': [0],
        'zone_count': [1],
        'fragmentation_count': [np.nan],
        'capture_count': [np.nan],
        'area_captured_sum': [np.nan],
        'area_captured_max': [np.nan]}
    )
    pd.testing.assert_frame_equal(
        sc.record_data_frame, expected_df, check_like=True
    )

    np.testing.assert_equal(len(se.filter_species(time=0)), 1)

    # Break the zone in two for time 1.

    z[11] = 0

    sc.run_one_step(1)
    se.run_one_step(1)

    np.testing.assert_equal(len(sc.zones), 2)
    np.testing.assert_equal(
        set([z._conn_type for z in sc.zones]), set([None, zn._ONE_TO_MANY])
    )

    expected_df = pd.DataFrame({
        'time': [0, 1],
        'zone_count': [1, 2],
        'fragmentation_count': [np.nan, 2],
        'capture_count': [np.nan, 0],
        'area_captured_sum': [np.nan, 0],
        'area_captured_max': [np.nan, 0]}
    )
    pd.testing.assert_frame_equal(
        sc.record_data_frame, expected_df, check_like=True
    )

    np.testing.assert_equal(len(se.filter_species(time=1)), 2)


def test_many_to_one(zone_example_grid):
    mg, z = zone_example_grid

    # Create two zones for time 0.

    z[[2, 3, 4, 17, 18]] = 1

    se = SpeciesEvolver(mg)
    sc = ZoneController(mg, zone_func)
    species = sc.populate_zones_uniformly(1)
    se.introduce_species(species)

    expected_df = pd.DataFrame({
        'time': [0],
        'zone_count': [2],
        'fragmentation_count': [np.nan],
        'capture_count': [np.nan],
        'area_captured_sum': [np.nan],
        'area_captured_max': [np.nan]}
    )
    pd.testing.assert_frame_equal(
        sc.record_data_frame, expected_df, check_like=True
    )

    np.testing.assert_equal(len(se.filter_species(time=0)), 2)

    # Modify elevation such that two zones each overlap the original two zones.

    z[[2, 3, 17]] = 0
    z[[11, 18]] = 1

    sc.run_one_step(1)
    se.run_one_step(1)

    np.testing.assert_equal(len(sc.zones), 1)
    np.testing.assert_equal(sc.zones[0]._conn_type, zn._MANY_TO_ONE)

    expected_df = pd.DataFrame({
        'time': [0, 1],
        'zone_count': [2, 1],
        'fragmentation_count': [np.nan, 0],
        'capture_count': [np.nan, 1],
        'area_captured_sum': [np.nan, 12],
        'area_captured_max': [np.nan, 12]}
    )
    pd.testing.assert_frame_equal(
        sc.record_data_frame, expected_df, check_like=True
    )

    np.testing.assert_equal(len(se.filter_species(time=1)), 2)


def test_many_to_many(zone_example_grid):
    mg, z = zone_example_grid

    # Create two zones for time 0.

    z[[2, 3, 4, 16, 17, 18]] = 1

    se = SpeciesEvolver(mg)
    sc = ZoneController(mg, zone_func)
    species = sc.populate_zones_uniformly(1)
    se.introduce_species(species)

    expected_df = pd.DataFrame({
        'time': [0],
        'zone_count': [2],
        'fragmentation_count': [np.nan],
        'capture_count': [np.nan],
        'area_captured_sum': [np.nan],
        'area_captured_max': [np.nan]}
    )
    pd.testing.assert_frame_equal(
        sc.record_data_frame, expected_df, check_like=True
    )

    np.testing.assert_equal(len(se.filter_species(time=0)), 2)

    # Modify elevation such that two zones each overlap the original two zones.

    z[[3, 17]] = 0
    z[[9, 11]] = 1

    sc.run_one_step(1)
    se.run_one_step(1)

    np.testing.assert_equal(len(sc.zones), 2)
    for z in sc.zones:
        np.testing.assert_equal(z._conn_type, zn._MANY_TO_MANY)

    expected_df = pd.DataFrame({
        'time': [0, 1],
        'zone_count': [2, 2],
        'fragmentation_count': [np.nan, 0],
        'capture_count': [np.nan, 2],
        'area_captured_sum': [np.nan, 24],
        'area_captured_max': [np.nan, 12]}
    )
    pd.testing.assert_frame_equal(
        sc.record_data_frame, expected_df, check_like=True
    )

    np.testing.assert_equal(len(se.filter_species(time=1)), 4)


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

    sc = ZoneController(mg, zone_func, neighborhood_structure='D4')
    np.testing.assert_equal(len(sc.zones), 2)

    np.testing.assert_raises(
        ValueError, ZoneController, mg, zone_func, neighborhood_structure='D'
    )


def test_zone_func_kwargs(zone_example_grid):
    mg, z = zone_example_grid
    z[:] = 1
    sc = ZoneController(mg, zone_func_with_vars, var1=1, var2=2)

    expected_mask = np.array(
        [False, False, False, False, False, False, False, False, False, False,
         False, False, False, False, False, True, True, True, True, True, True]
    )
    np.testing.assert_array_equal(sc.zones[0].mask, expected_mask)

    z[-2:] = 0
    expected_mask[-2:] = False

    sc.run_one_step(10)

    np.testing.assert_array_equal(sc.zones[0].mask, expected_mask)


def test_zone_species_range_mask(zone_example_grid):
    mg, z = zone_example_grid

    # Create a zone for time 0.

    ids = [9, 10, 11, 12]
    z[ids] = 1

    se = SpeciesEvolver(mg)
    sc = ZoneController(mg, zone_func)
    species = sc.populate_zones_uniformly(1)
    se.introduce_species(species)

    expected_mask = np.zeros(mg.number_of_nodes, bool)
    expected_mask[ids] = True
    np.testing.assert_array_equal(species[0].range_mask, expected_mask)

    # Remove extent so species range mask is all False.

    z[ids] = 0
    sc.run_one_step(1)
    se.run_one_step(1)

    expected_mask = np.zeros(mg.number_of_nodes, bool)
    np.testing.assert_array_equal(species[0].range_mask, expected_mask)


def test_allopatric_wait_time(zone_example_grid):
    mg, z = zone_example_grid

    # Create a zone for time 0.

    z[[9, 10, 11, 12]] = 1

    se = SpeciesEvolver(mg)
    sc = ZoneController(mg, zone_func)
    species = sc.populate_zones_uniformly(1, allopatric_wait_time=20)
    se.introduce_species(species)

    z[[11]] = 0

    while len(se.species_data_frame) == 1:
        sc.run_one_step(10)
        se.run_one_step(10)

    expected_df = pd.DataFrame({
        'clade': ['A', 'A', 'A'],
        'number': [0, 1, 2],
        'time_appeared': [0, 30, 30],
        'latest_time': [20, 30, 30]}
    )
    pd.testing.assert_frame_equal(
        se.species_data_frame, expected_df, check_like=True
    )

#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""SpeciesEvolver tests."""
from collections import Counter

import numpy as np
import pandas as pd
import pytest

from landlab import RasterModelGrid
from landlab.components import SpeciesEvolver
from landlab.components.species_evolution.base_species import Species


class SpeciesTest(Species):
    """Species to test SpeciesEvolver."""

    def __init__(self, parent_species=None):
        super(SpeciesTest, self).__init__()
        self._parent_species = parent_species

    @property
    def range_mask(self):
        return np.array(
            [False, False, False, True, True, True, False, False, False]
        )

    def _evolve_stage_1(self, dt, record):
        pass

    def _evolve_stage_2(self, dt, record):
        species_persists = False
        child_species = [SpeciesTest(parent_species=self)]
        return species_persists, child_species


@pytest.fixture()
def zone_example_grid():
    return RasterModelGrid((3, 3), 1)


def test_introduce_species_and_component_attributes(zone_example_grid):
    se = SpeciesEvolver(zone_example_grid)

    # Introduce multiple species.
    species = [SpeciesTest(), SpeciesTest()]
    se.introduce_species(species)

    # Introduce a single species.
    species = SpeciesTest()
    se.introduce_species(species)

    # Test attributes at initial time step.

    expected_df = pd.DataFrame({
        'clade': ['A', 'B', 'C'],
        'number': [0, 0, 0],
        'time_appeared': [0, 0, 0],
        'latest_time': [0, 0, 0]
    })
    pd.testing.assert_frame_equal(
        se.species_data_frame, expected_df, check_like=True
    )

    expected_df = pd.DataFrame({
        'time': [0],
        'species_count': [3]
    })
    pd.testing.assert_frame_equal(
        se.record_data_frame, expected_df, check_like=True
    )

    # Test attributes at a later time.

    se.run_one_step(10)

    expected_df = pd.DataFrame({
        'clade': ['A', 'A', 'B', 'B', 'C', 'C'],
        'number': [0, 1, 0, 1, 0, 1],
        'time_appeared': [0, 10, 0, 10, 0, 10],
        'latest_time': [0, 10, 0, 10, 0, 10]
    })
    pd.testing.assert_frame_equal(
        se.species_data_frame, expected_df, check_like=True
    )

    expected_df = pd.DataFrame({
        'time': [0, 10],
        'species_count': [3, 3]
    })
    pd.testing.assert_frame_equal(
        se.record_data_frame, expected_df, check_like=True
    )


def test_filter_species(zone_example_grid):
    se = SpeciesEvolver(zone_example_grid)

    introduced_species = [SpeciesTest(), SpeciesTest()]
    se.introduce_species(introduced_species)
    se.run_one_step(10)

    # Test no parameters.

    queried_species = se.filter_species()
    np.testing.assert_equal(
        Counter(queried_species), Counter(se._species['object'])
    )

    # Test `time` parameter.

    queried_species = se.filter_species(time=0)
    np.testing.assert_equal(
        Counter(queried_species), Counter(introduced_species)
    )

    queried_species = se.filter_species(time=10)
    child_species = set(se._species['object']) - set(introduced_species)
    np.testing.assert_equal(
        Counter(queried_species), Counter(child_species)
    )

    np.testing.assert_raises(ValueError, se.filter_species, time=5)
    np.testing.assert_raises(ValueError, se.filter_species, time=11)

    # Test `clade` parameter.

    queried_species = se.filter_species(clade='B')
    ids = [s.identifier for s in queried_species]
    expected_ids = [('B', 0), ('B', 1)]
    np.testing.assert_equal(Counter(ids), Counter(expected_ids))

    queried_species = se.filter_species(clade='C')
    np.testing.assert_equal(queried_species, [])

    # Test `number` parameter.

    queried_species = se.filter_species(number=1)
    np.testing.assert_equal(
        Counter(queried_species), Counter(child_species)
    )

    queried_species = se.filter_species(time=10, number=3)
    np.testing.assert_equal(queried_species, [])

    # Test multiple parameters.

    queried_species = se.filter_species(clade='B', number=1)
    np.testing.assert_equal(
        queried_species[0].identifier, ('B', 1)
    )

    queried_species = se.filter_species(time=10, number=1)
    ids = [s.identifier for s in queried_species]
    expected_ids = [('A', 1), ('B', 1)]
    np.testing.assert_equal(Counter(ids), Counter(expected_ids))

    queried_species = se.filter_species(time=10, number=0)
    np.testing.assert_equal(queried_species, [])

    queried_species = se.filter_species(time=10, clade='B', number=1)
    np.testing.assert_equal(
        queried_species[0].identifier, ('B', 1)
    )

    # Test `species_subset` parameter.

    queried_species_1 = se.filter_species(time=0)
    queried_species_2 = se.filter_species(species_subset=queried_species_1)
    np.testing.assert_equal(
        Counter(queried_species_1), Counter(queried_species_2)
    )
    queried_species_3 = se.filter_species(
        species_subset=queried_species_1, time=10
    )
    np.testing.assert_equal(queried_species_3, [])
    queried_species_4 = se.filter_species(
        species_subset=queried_species_1, clade='B'
    )
    np.testing.assert_equal(queried_species_4[0].identifier, ('B', 0))


def test_species_richness_field(zone_example_grid):
    mg = zone_example_grid

    se = SpeciesEvolver(mg)

    expected_field = np.zeros(mg.number_of_nodes)
    np.testing.assert_array_equal(
        mg.at_node['species__richness'], expected_field
    )

    introduced_species = [SpeciesTest(), SpeciesTest()]
    se.introduce_species(introduced_species)
    se.run_one_step(10)

    expected_field = np.array([0, 0, 0, 2, 2, 2, 0, 0, 0])
    np.testing.assert_array_equal(
        mg.at_node['species__richness'], expected_field
    )

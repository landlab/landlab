#!/usr/bin/env python
"""SpeciesEvolver tests."""
from collections import Counter

import numpy as np
from pandas import DataFrame
import pytest

from landlab import RasterModelGrid
from landlab.components import SpeciesEvolver
from landlab.components.species_evolution import _Species


class SpeciesTest(_Species):
    """Species to test SpeciesEvolver."""

    def __init__(self, parent_species=None):
        super(SpeciesTest, self).__init__()
        self._parent_species = parent_species

    @property
    def range_mask(self):
        pass

    def _evolve_stage_1(self, dt, record):
        pass

    def _evolve_stage_2(self, dt, record):
        species_persists = True
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

    expected_df = DataFrame({
        'clade': ['A', 'B', 'C'],
        'number': [0, 0, 0],
        'time_appeared': [0, 0, 0],
        'latest_time': [0, 0, 0]
    })
    np.testing.assert_array_equal(se.species_data_frame, expected_df)

    expected_df = DataFrame({
        'time': [0],
        'species_count': [3]
    })
    np.testing.assert_array_equal(se.record_data_frame, expected_df)

    # Test attributes at a later time.

    se.run_one_step(10)

    expected_df = DataFrame({
        'clade': ['A', 'A', 'B', 'B', 'C', 'C'],
        'number': [0, 1, 0, 1, 0, 1],
        'time_appeared': [0, 10, 0, 10, 0, 10],
        'latest_time': [10, 10, 10, 10, 10, 10]
    })
    np.testing.assert_array_equal(se.species_data_frame, expected_df)

    expected_df = DataFrame({
        'time': [0, 10],
        'species_count': [3, 6]
    })
    np.testing.assert_array_equal(se.record_data_frame, expected_df)


def test_species_at_time(zone_example_grid):
    se = SpeciesEvolver(zone_example_grid)
    introduced_species = [SpeciesTest(), SpeciesTest()]
    se.introduce_species(introduced_species)
    se.run_one_step(10)

    # Test time steps in the SpeciesEvolver record.

    queried_species = se.species_at_time(0)
    np.testing.assert_equal(Counter(queried_species),
        Counter(introduced_species))

    queried_species = se.species_at_time(10)
    ids = [s.identifier for s in queried_species]
    expected_ids = [('A', 0), ('A', 1), ('B', 0), ('B', 1)]
    np.testing.assert_equal(Counter(ids), Counter(expected_ids))

    # Test time steps in between and outside of the SpeciesEvolver record.

    queried_species = se.species_at_time(5)
    np.testing.assert_equal(Counter(queried_species),
        Counter(introduced_species))

    queried_species = se.species_at_time(-1)
    np.testing.assert_equal(queried_species, np.nan)

    queried_species = se.species_at_time(11)
    np.testing.assert_equal(queried_species, np.nan)


def test_species_with_identifier(zone_example_grid):
    se = SpeciesEvolver(zone_example_grid)
    introduced_species = [SpeciesTest(), SpeciesTest()]
    se.introduce_species(introduced_species)
    se.run_one_step(10)

    queried_species = se.species_with_identifier(('B', 1))
    np.testing.assert_equal(queried_species[0].identifier, ('B', 1))

    queried_species = se.species_with_identifier('B')
    ids = [s.identifier for s in queried_species]
    expected_ids = [('B', 0), ('B', 1)]
    np.testing.assert_equal(Counter(ids), Counter(expected_ids))

    queried_species = se.species_with_identifier(1)
    ids = [s.identifier for s in queried_species]
    expected_ids = [('A', 1), ('B', 1)]
    np.testing.assert_equal(Counter(ids), Counter(expected_ids))

    np.testing.assert_raises(TypeError, se.species_with_identifier, (1, 'A'))

    np.testing.assert_raises(TypeError, se.species_with_identifier, {})

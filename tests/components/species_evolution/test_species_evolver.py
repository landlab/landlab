#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""SpeciesEvolver tests."""
from collections import Counter

import numpy as np
import pandas as pd
import pytest

from landlab import RasterModelGrid
from landlab.components import SpeciesEvolver
from landlab.components.species_evolution.base_taxon import Taxon


class TestTaxon(Taxon):
    """Taxon to test SpeciesEvolver."""

    def __init__(self, parent=None):
        super(TestTaxon, self).__init__()
        self.parent = parent

    @property
    def range_mask(self):
        return np.array(
            [False, False, False, True, True, True, False, False, False]
        )

    def _evolve(self, dt, stage, record):
        if stage == 1:
            self._extant = False
            TestTaxon(parent=self)

        return stage < 1


@pytest.fixture()
def zone_example_grid():
    return RasterModelGrid((3, 3), 1)


def test_track_taxa_and_component_attributes(zone_example_grid):
    se = SpeciesEvolver(zone_example_grid)

    # Introduce multiple taxa.
    taxa = [TestTaxon(), TestTaxon()]
    se.track_taxa(taxa)

    # Introduce a single taxon.
    taxon = TestTaxon()
    se.track_taxa(taxon)

    # Test attributes at initial time step.

    expected_df = pd.DataFrame({
        'appeared': [0, 0, 0],
        'latest_time': [0, 0, 0],
        'extant': [True, True, True]},
        index=[0, 1, 2]
    )
    expected_df.index.name = 'uid'
    pd.testing.assert_frame_equal(
        se.taxa_data_frame, expected_df, check_like=True
    )

    expected_df = pd.DataFrame({
        'time': [0],
        'taxa': [3]
    })
    pd.testing.assert_frame_equal(
        se.record_data_frame, expected_df, check_like=True
    )

    # Test attributes at a later time.

    se.run_one_step(10)

    expected_df = pd.DataFrame({
        'appeared': [0, 0, 0, 10, 10, 10],
        'latest_time': [10, 10, 10, 10, 10, 10],
        'extant': [False, False, False, True, True, True]},
        index=[0, 1, 2, 3, 4, 5]
    )
    pd.testing.assert_frame_equal(
        se.taxa_data_frame, expected_df, check_like=True
    )

    expected_df = pd.DataFrame({
        'time': [0, 10],
        'taxa': [3, 3]
    })
    pd.testing.assert_frame_equal(
        se.record_data_frame, expected_df, check_like=True
    )


def test_get_taxon_objects(zone_example_grid):
    se = SpeciesEvolver(zone_example_grid)

    introduced_taxa = [TestTaxon(), TestTaxon()]
    se.track_taxa(introduced_taxa)
    se.run_one_step(10)
    se.run_one_step(10)

    # Test no parameters.

    queried_taxa = se.get_taxon_objects()
    np.testing.assert_equal(
        Counter(queried_taxa), Counter(se._taxa['object'])
    )

    # Test `time` parameter.

    queried_taxa = se.get_taxon_objects(time=0)
    np.testing.assert_equal(
        Counter(queried_taxa), Counter(introduced_taxa)
    )

    queried_taxa = se.get_taxon_objects(time=10)
    ids = [s.uid for s in queried_taxa]
    expected_ids = [0, 1, 2, 3]
    np.testing.assert_equal(Counter(ids), Counter(expected_ids))

    np.testing.assert_raises(ValueError, se.get_taxon_objects, time=5)
    np.testing.assert_raises(ValueError, se.get_taxon_objects, time=11)

    # Test `extant_at_latest_time` parameter.

    queried_taxa = se.get_taxon_objects(extant_at_latest_time=True)
    ids = [s.uid for s in queried_taxa]
    expected_ids = [4, 5]
    np.testing.assert_equal(Counter(ids), Counter(expected_ids))

    # Test `ancestor` parameter.

    queried_taxa = se.get_taxon_objects(ancestor=1)
    ids = [s.uid for s in queried_taxa]
    expected_ids = [3, 5]
    np.testing.assert_equal(Counter(ids), Counter(expected_ids))

    queried_taxa = se.get_taxon_objects(ancestor=5)
    np.testing.assert_equal(queried_taxa, [])

    queried_taxa = se.get_taxon_objects(ancestor=6)
    np.testing.assert_equal(queried_taxa, [])

    # Test multiple parameters.

    queried_taxa = se.get_taxon_objects(ancestor=1, time=10)
    ids = [s.uid for s in queried_taxa]
    expected_ids = [3]
    np.testing.assert_equal(Counter(ids), Counter(expected_ids))


def test_taxa_richness_field(zone_example_grid):
    mg = zone_example_grid

    se = SpeciesEvolver(mg)

    expected_field = np.zeros(mg.number_of_nodes)
    np.testing.assert_array_equal(
        mg.at_node['taxa__richness'], expected_field
    )

    introduced_taxa = [TestTaxon(), TestTaxon()]
    se.track_taxa(introduced_taxa)
    se.run_one_step(10)

    expected_field = np.array([0, 0, 0, 2, 2, 2, 0, 0, 0])
    np.testing.assert_array_equal(
        mg.at_node['taxa__richness'], expected_field
    )

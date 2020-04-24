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
        super().__init__()
        self.parent = parent

    @property
    def range_mask(self):
        return np.array([False, False, False, True, True, True, False, False, False])

    def _evolve(self, dt, stage, record):
        if stage == 0:
            self._extant = False
            child_taxon = [TestTaxon(parent=self)]

        return False, child_taxon


@pytest.fixture()
def zone_example_grid():
    return RasterModelGrid((3, 3))


def test_track_taxa_and_component_attributes(zone_example_grid):
    se = SpeciesEvolver(zone_example_grid)

    # Introduce multiple taxa.
    taxa = [TestTaxon(), TestTaxon()]
    se.track_taxa(taxa)

    # Introduce a single taxon.
    taxon = TestTaxon()
    se.track_taxa(taxon)

    # Test attributes at initial time step.

    expected_df = pd.DataFrame(
        {
            "pid": 3 * [np.nan],
            "type": 3 * [TestTaxon.__name__],
            "t_first": [0, 0, 0],
            "t_final": 3 * [np.nan],
        },
        index=[0, 1, 2],
    )
    expected_df.index.name = "tid"
    expected_df["pid"] = expected_df["pid"].astype("Int64")
    expected_df["t_final"] = expected_df["t_final"].astype("Int64")
    pd.testing.assert_frame_equal(se.taxa_data_frame, expected_df, check_like=True)

    expected_df = pd.DataFrame({"time": [0], "taxa": [3]})
    pd.testing.assert_frame_equal(se.record_data_frame, expected_df, check_like=True)

    # Test attributes at a later time.

    se.run_one_step(10)

    expected_df = pd.DataFrame(
        {
            "pid": 3 * [np.nan] + [0, 1, 2],
            "type": 6 * [TestTaxon.__name__],
            "t_first": [0, 0, 0, 10, 10, 10],
            "t_final": [10, 10, 10] + 3 * [np.nan],
        },
        index=[0, 1, 2, 3, 4, 5],
    )
    expected_df.index.name = "tid"
    expected_df["pid"] = expected_df["pid"].astype("Int64")
    expected_df["t_final"] = expected_df["t_final"].astype("Int64")
    pd.testing.assert_frame_equal(se.taxa_data_frame, expected_df, check_like=True)

    expected_df = pd.DataFrame({"time": [0, 10], "taxa": [3, 3]})
    pd.testing.assert_frame_equal(se.record_data_frame, expected_df, check_like=True)


def test_get_extant_taxon_objects(zone_example_grid):
    se = SpeciesEvolver(zone_example_grid)

    introduced_taxa = [TestTaxon(), TestTaxon()]
    se.track_taxa(introduced_taxa)
    se.run_one_step(10)
    se.run_one_step(10)

    # Test no parameters.

    queried_taxa = se.get_extant_taxon_objects()
    np.testing.assert_equal(Counter(queried_taxa), Counter(se._taxon_objs))

    # Test `tids` parameter.

    queried_taxa = se.get_extant_taxon_objects(tids=[0])
    ids = [t.tid for t in queried_taxa]
    expected_ids = []
    np.testing.assert_equal(Counter(ids), Counter(expected_ids))

    queried_taxa = se.get_extant_taxon_objects(tids=[4, 5])
    ids = [t.tid for t in queried_taxa]
    expected_ids = [4, 5]
    np.testing.assert_equal(Counter(ids), Counter(expected_ids))

    # Test `time` parameter.

    queried_taxa = se.get_extant_taxon_objects(time=20)
    ids = [t.tid for t in queried_taxa]
    expected_ids = [4, 5]
    np.testing.assert_equal(Counter(ids), Counter(expected_ids))

    queried_taxa = se.get_extant_taxon_objects(time=10)
    ids = [t.tid for t in queried_taxa]
    expected_ids = []
    np.testing.assert_equal(Counter(ids), Counter(expected_ids))

    queried_taxa = se.get_extant_taxon_objects(time=30)
    ids = [t.tid for t in queried_taxa]
    expected_ids = []
    np.testing.assert_equal(Counter(ids), Counter(expected_ids))

    # Test `ancestor` parameter.

    queried_taxa = se.get_extant_taxon_objects(ancestor=1)
    ids = [t.tid for t in queried_taxa]
    expected_ids = [5]
    np.testing.assert_equal(Counter(ids), Counter(expected_ids))

    queried_taxa = se.get_extant_taxon_objects(ancestor=5)
    np.testing.assert_equal(queried_taxa, [])

    queried_taxa = se.get_extant_taxon_objects(ancestor=6)
    np.testing.assert_equal(queried_taxa, [])

    # Test multiple parameters.

    queried_taxa = se.get_extant_taxon_objects(ancestor=1, time=10)
    np.testing.assert_equal(queried_taxa, [])


def test_taxa_richness_field(zone_example_grid):
    mg = zone_example_grid

    se = SpeciesEvolver(mg)

    expected_field = np.zeros(mg.number_of_nodes)
    np.testing.assert_array_equal(mg.at_node["taxa__richness"], expected_field)

    introduced_taxa = [TestTaxon(), TestTaxon()]
    se.track_taxa(introduced_taxa)
    se.run_one_step(10)

    expected_field = np.array([0, 0, 0, 2, 2, 2, 0, 0, 0])
    np.testing.assert_array_equal(mg.at_node["taxa__richness"], expected_field)


def test_immediate_extinction(zone_example_grid):
    se = SpeciesEvolver(zone_example_grid)

    taxon = TestTaxon()
    taxon.extant = False
    se.track_taxa(taxon)

    expected_df = pd.DataFrame(
        {"pid": [np.nan], "type": [TestTaxon.__name__], "t_first": [0], "t_final": [0]},
        index=[0],
    )
    expected_df.index.name = "tid"
    expected_df["pid"] = expected_df["pid"].astype("Int64")
    expected_df["t_final"] = expected_df["t_final"].astype("Int64")
    pd.testing.assert_frame_equal(se.taxa_data_frame, expected_df, check_like=True)

#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""Tests for the base species of SpeciesEvolver."""
from collections import OrderedDict

import numpy as np

from landlab.components.species_evolution.base_species import Species
from landlab.components.species_evolution.record import Record


class SpeciesTest(Species):
    """Species to test SpeciesEvolver."""

    def __init__(self, parent_species=None):
        super(SpeciesTest, self).__init__()
        self._parent_species = parent_species

    @property
    def range_mask(self):
        return True

    def _evolve_stage_1(self, dt, record):
        record.set_value('vara', 1)

    def _evolve_stage_2(self, dt, record):
        species_persists = True
        child_species = [SpeciesTest(parent_species=self)]
        record.set_value('varb', 2)
        return species_persists, child_species


def test_base():
    ts = SpeciesTest()
    assert isinstance(ts, Species)
    np.testing.assert_equal(ts.range_mask, True)

    record = Record()
    record.advance_time(10)

    output = ts._evolve_stage_1(10, record)
    np.testing.assert_equal(output, None)

    d = OrderedDict([('time', [0, 10]), ('vara', [np.nan, 1])])
    np.testing.assert_equal(record._dict, d)

    species_persists, child_species = ts._evolve_stage_2(10, record)
    np.testing.assert_equal(species_persists, True)
    np.testing.assert_equal(type(child_species[0]), SpeciesTest)

    d = OrderedDict(
        [('time', [0, 10]), ('vara', [np.nan, 1]), ('varb', [np.nan, 2])]
    )
    np.testing.assert_equal(record._dict, d)

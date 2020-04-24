#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""Tests for the base taxon of SpeciesEvolver."""
from collections import OrderedDict

import numpy as np

from landlab.components.species_evolution.base_taxon import Taxon
from landlab.components.species_evolution.record import Record


class TestTaxon(Taxon):
    """Taxon to test subclassing the base taxon."""

    def __init__(self, parent=None):
        super().__init__()
        self.parent = parent

    @property
    def range_mask(self):
        return True

    def _evolve(self, dt, stage, record):
        if stage == 0:
            record.set_value("vara", 1)
        elif stage == 1:
            TestTaxon(parent=self)
            record.set_value("varb", 2)

        return stage < 1


def test_base_subclass():
    tt = TestTaxon()
    assert isinstance(tt, Taxon)
    assert tt.__repr__() == "<TestTaxon, tid=None>"

    np.testing.assert_equal(tt.range_mask, True)

    record = Record()
    record.advance_time(10)

    output = tt._evolve(10, 0, record)
    np.testing.assert_equal(output, True)

    d = OrderedDict([("time", [0, 10]), ("vara", [np.nan, 1])])
    np.testing.assert_equal(record._dict, d)

    output = tt._evolve(10, 1, record)
    np.testing.assert_equal(output, False)

    np.testing.assert_equal(tt.extant, True)

    d = OrderedDict([("time", [0, 10]), ("vara", [np.nan, 1]), ("varb", [np.nan, 2])])
    np.testing.assert_equal(record._dict, d)

    tt.extant = False
    np.testing.assert_equal(tt._extant, False)

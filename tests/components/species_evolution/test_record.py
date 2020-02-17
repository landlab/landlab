#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""Tests for Record of SpeciesEvolver."""
from collections import OrderedDict

import numpy as np
import pandas as pd

from landlab.components.species_evolution.record import Record


def test_attributes():
    record = Record()
    np.testing.assert_equal(record.prior_time, np.nan)

    record.advance_time(10)

    record.set_value("vara", 1)
    record.set_value("varb", 2)

    np.testing.assert_equal(record.times, [0, 10])
    np.testing.assert_equal(record.prior_time, 0)
    np.testing.assert_equal(record.earliest_time, 0)
    np.testing.assert_equal(record.latest_time, 10)
    np.testing.assert_equal(record.variables, ["vara", "varb"])
    np.testing.assert_equal(record.count_of_time_steps, 2)
    np.testing.assert_equal(len(record), 2)

    df = pd.DataFrame({"time": [0, 10], "vara": [np.nan, 1], "varb": [np.nan, 2]})
    pd.testing.assert_frame_equal(record.data_frame, df, check_like=True)


def test_get_value():
    record = Record()
    record.set_value("vara", 1)
    record.advance_time(10)
    record.set_value("vara", 2)

    value = record.get_value("vara")
    np.testing.assert_equal(value, 2)

    value = record.get_value("vara", time=0)
    np.testing.assert_equal(value, 1)

    value = record.get_value("varb")
    np.testing.assert_equal(value, np.nan)


def test_set_value():
    record = Record()
    record.advance_time(10)

    record.set_value("vara", 1)
    d = OrderedDict([("time", [0, 10]), ("vara", [np.nan, 1])])
    np.testing.assert_equal(record._dict, d)

    record.set_value("vara", 3, time=0)
    d = OrderedDict([("time", [0, 10]), ("vara", [3, 1])])
    np.testing.assert_equal(record._dict, d)


def test_increment_value():
    record = Record()
    record.set_value("vara", 1)
    record.increment_value("vara", 100)

    d = OrderedDict([("time", [0]), ("vara", [101])])
    np.testing.assert_equal(record._dict, d)

    record.advance_time(10)
    record.increment_value("vara", 100)

    d = OrderedDict([("time", [0, 10]), ("vara", [101, 100])])
    np.testing.assert_equal(record._dict, d)

    record.increment_value("varb", 200, time=0)

    d = OrderedDict([("time", [0, 10]), ("vara", [101, 100]), ("varb", [200, np.nan])])
    np.testing.assert_equal(record._dict, d)


def test_initial_time():
    record = Record(initial_time=5)
    record.set_value("vara", 1)
    record.advance_time(10)

    d = OrderedDict([("time", [5, 15]), ("vara", [1, np.nan])])
    np.testing.assert_equal(record._dict, d)


def test_nonexistent_time():
    record = Record()
    record.set_value("vara", 1)

    np.testing.assert_raises(ValueError, record.get_value, "vara", 5)

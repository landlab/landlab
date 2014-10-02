#! /usr/bin/env python

from nose.tools import assert_equal
import numpy as np

from landlab.utils import count_repeated_values


def test_empty_array():
    x = np.array([])
    counts = count_repeated_values(x)
    assert_equal(len(counts), 0)


def test_no_dups():
    x = np.array([10, 20, 30])
    counts = count_repeated_values(x)
    assert_equal(len(counts), 1)

    (vals, inds) = counts[0]
    assert_equal(list(vals), [10, 20, 30])
    assert_equal(list(inds), [0, 1, 2])


def test_with_dups():
    x = np.array([10, 20, 30, 20])
    counts = count_repeated_values(x)
    assert_equal(len(counts), 2)

    (vals, inds) = counts[0]
    assert_equal(list(vals), [10, 20, 30])
    assert_equal(list(inds), [0, 1, 2])

    (vals, inds) = counts[1]
    assert_equal(list(vals), [20, ])
    assert_equal(list(inds), [3, ])

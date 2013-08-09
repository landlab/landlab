#! /usr/bin/env python

import unittest
import numpy as np

from landlab.utils import count_repeated_values


class TestCountRepeats(unittest.TestCase):
    def test_empty_array(self):
        x = np.array([])
        counts = count_repeated_values(x)
        self.assertEqual(len(counts), 0)

    def test_no_dups(self):
        x = np.array([10, 20, 30])
        counts = count_repeated_values(x)
        self.assertEqual(len(counts), 1)

        (vals, inds) = counts[0]
        self.assertEqual(list(vals), [10, 20, 30])
        self.assertEqual(list(inds), [0, 1, 2])

    def test_with_dups(self):
        x = np.array([10, 20, 30, 20])
        counts = count_repeated_values(x)
        self.assertEqual(len(counts), 2)

        (vals, inds) = counts[0]
        self.assertEqual(list(vals), [10, 20, 30])
        self.assertEqual(list(inds), [0, 1, 2])

        (vals, inds) = counts[1]
        self.assertEqual(list(vals), [20, ])
        self.assertEqual(list(inds), [3, ])


if __name__ == '__main__':
    unittest.main()

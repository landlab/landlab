#! /usr/bin/env python

import unittest
import numpy as np
from numpy.testing import assert_array_equal

from landlab.field.grouped import ModelDataFields


class TestModelDataFields(unittest.TestCase):
    def test_init(self):
        fields = ModelDataFields()
        self.assertSetEqual(set(), fields.groups)

    def test_new_field_location(self):
        fields = ModelDataFields()
        fields.new_field_location('node', 12)
        self.assertSetEqual(set(['node']), fields.groups)

    def test_add_existing_group(self):
        fields = ModelDataFields()
        fields.new_field_location('node', 12)
        with self.assertRaises(ValueError):
            fields.new_field_location('node', 24)

    def test_add_multiple_groups(self):
        fields = ModelDataFields()
        fields.new_field_location('node', 12)
        fields.new_field_location('cell', 2)
        fields.new_field_location('face', 7)
        fields.new_field_location('link', 7)
        self.assertSetEqual(set(['node', 'cell', 'face', 'link']),
                            fields.groups)

    def test_ones(self):
        fields = ModelDataFields()
        fields.new_field_location('node', 12)
        fields.new_field_location('cell', 2)

        value_array = fields.ones('node')
        assert_array_equal(np.ones(12), value_array)

        value_array = fields.ones('cell')
        assert_array_equal(np.ones(2), value_array)

    def test_add_ones(self):
        fields = ModelDataFields()
        fields.new_field_location('node', 12)
        fields.new_field_location('cell', 2)

        fields.add_ones('node', 'z')
        assert_array_equal(np.ones(12), fields['node']['z'])
        assert_array_equal(np.ones(12), fields.field_values('node', 'z'))

        fields.add_ones('cell', 'z')
        assert_array_equal(np.ones(2), fields['cell']['z'])
        assert_array_equal(np.ones(2), fields.field_values('cell', 'z'))

    def test_add_ones_return_value(self):
        fields = ModelDataFields()
        fields.new_field_location('node', 12)
        fields.new_field_location('cell', 2)

        rtn_value = fields.add_ones('node', 'z')
        assert_array_equal(rtn_value, np.ones(12))
        self.assertIs(rtn_value, fields['node']['z'])
        self.assertIs(rtn_value, fields.field_values('node', 'z'))

        rtn_value = fields.add_ones('cell', 'z')
        assert_array_equal(rtn_value, np.ones(2))
        self.assertIs(rtn_value, fields['cell']['z'])
        self.assertIs(rtn_value, fields.field_values('cell', 'z'))

    def test_getitem(self):
        fields = ModelDataFields()
        fields.new_field_location('node', 12)

        self.assertDictEqual(dict(), fields['node'])
        with self.assertRaises(KeyError):
            fields['cell']

    def test_at_attribute(self):
        fields = ModelDataFields()
        fields.new_field_location('node', 12)

        self.assertDictEqual(dict(), fields.at_node)
        with self.assertRaises(AttributeError):
            fields.at_cell

        fields.add_ones('node', 'z')
        assert_array_equal(np.ones(12), fields.at_node['z'])

    def test_has_group(self):
        fields = ModelDataFields()
        fields.new_field_location('node', 12)

        self.assertTrue(fields.has_group('node'))
        self.assertFalse(fields.has_group('cell'))

#! /usr/bin/env python
"""
Unit tests for landlab.model_parameter_dictionary
"""

import unittest

from landlab import ModelParameterDictionary
from landlab.model_parameter_dictionary import (MissingKeyError,
                                                ParameterValueError)

_TEST_PARAM_DICT_FILE = """
# A Comment
INT_VAL:
1

FLOAT_VAL: # A Comment
2.2
STRING_VAL:
The Landlab

TRUE_BOOL_VAL:
True
FALSE_BOOL_VAL:
False
"""

class TestModelParameterDictionary(unittest.TestCase):
    def setUp(self):
        from StringIO import StringIO
        self.param_file = StringIO(_TEST_PARAM_DICT_FILE)
        self.param_dict = ModelParameterDictionary()

        self.param_dict.read_from_file(self.param_file)

    def test_read_file(self):
        all_keys = [
            'FLOAT_VAL', 'INT_VAL', 'STRING_VAL', 'TRUE_BOOL_VAL',
            'FALSE_BOOL_VAL',
        ]
        all_keys.sort()
        param_list = self.param_dict.params()
        param_list.sort()
        self.assertEqual(param_list, all_keys)

    def test_read_int(self):
        self.assertEqual(self.param_dict.read_int('INT_VAL'), 1)

        with self.assertRaises(ParameterValueError):
            self.param_dict.read_int('FLOAT_VAL')

        with self.assertRaises(MissingKeyError):
            self.param_dict.read_int('MISSING_INT')


    def test_read_float(self):
        self.assertEqual(self.param_dict.read_float('FLOAT_VAL'), 2.2)
        self.assertEqual(self.param_dict.read_float('INT_VAL'), 1)

        with self.assertRaises(ParameterValueError):
            self.param_dict.read_float('STRING_VAL')

        with self.assertRaises(MissingKeyError):
            self.param_dict.read_float('MISSING_FLOAT')

    def test_read_string(self):
        self.assertEqual(self.param_dict.read_string('STRING_VAL'),
                         'The Landlab')
        self.assertEqual(self.param_dict.read_string('INT_VAL'), '1')
        self.assertEqual(self.param_dict.read_string('FLOAT_VAL'), '2.2')

        with self.assertRaises(MissingKeyError):
            self.param_dict.read_string('MISSING_STRING')

    def test_read_bool(self):
        self.assertEqual(self.param_dict.read_bool('TRUE_BOOL_VAL'), True)
        self.assertEqual(self.param_dict.read_bool('FALSE_BOOL_VAL'), False)

        with self.assertRaises(MissingKeyError):
            self.param_dict.read_bool('MISSING_BOOLEAN')

        with self.assertRaises(ParameterValueError):
            self.param_dict.read_bool('STRING_VAL')


if __name__ == '__main__':
    unittest.main()

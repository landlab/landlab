#! /usr/bin/env python
"""
Unit tests for landlab.model_parameter_dictionary
"""

import unittest
import numpy as np

from landlab import ModelParameterDictionary
from landlab.model_parameter_dictionary import (MissingKeyError,
                                                ParameterValueError)

_TEST_PARAM_DICT_FILE = u"""
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
        from io import StringIO
        self.param_file = StringIO(_TEST_PARAM_DICT_FILE)
        self.param_dict = ModelParameterDictionary()

        self.param_dict.read_from_file(self.param_file)

    def test_read_file(self):
        all_keys = set([
            'FLOAT_VAL', 'INT_VAL', 'STRING_VAL', 'TRUE_BOOL_VAL',
            'FALSE_BOOL_VAL',
        ])
        param_list = set(self.param_dict.params())

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

    def test_dict_keys(self):
        all_keys = set([
            'FLOAT_VAL', 'INT_VAL', 'STRING_VAL', 'TRUE_BOOL_VAL',
            'FALSE_BOOL_VAL',
        ])
        self.assertEqual(set(self.param_dict), all_keys)
        for key in all_keys:
            self.assertTrue(key in self.param_dict)

    def test_dict_index(self):
        self.assertEqual(self.param_dict['INT_VAL'], '1')

class TestModelParameterDictionaryAutoType(unittest.TestCase):
    def setUp(self):
        from io import StringIO

        _TEST_FILE = u"""
# A Comment
INT_VAL:
1
DBL_VAL:
1.2
STR_VAL:
    landlab
INT_ARRAY_VAL:
1,2 ,4 ,7

DBL_ARRAY_VAL:
1.,2. ,4. ,7.
        """
        self.param_dict = ModelParameterDictionary(
            auto_type=True, from_file=StringIO(_TEST_FILE))

    def test_auto_type(self):
        self.assertEqual(self.param_dict['INT_VAL'], 1)
        self.assertEqual(self.param_dict['DBL_VAL'], 1.2)
        self.assertEqual(self.param_dict['STR_VAL'], 'landlab')

    def _test_int_vector(self):
        val = self.param_dict['INT_ARRAY_VAL']

        self.assertEqual(list(val), [1, 2, 4, 7])
        self.assertIsInstance(val, np.ndarray)
        self.assertEqual(val.dtype, np.int)

    def _test_float_vector(self):
        val = self.param_dict['DBL_ARRAY_VAL']

        self.assertEqual(list(val), [1., 2., 4., 7.])
        self.assertIsInstance(val, np.ndarray)
        self.assertEqual(val.dtype, np.float)


if __name__ == '__main__':
    unittest.main()

#! /usr/bin/env python
"""
Unit tests for landlab.model_parameter_dictionary
"""

from nose.tools import assert_equal, assert_true, assert_raises, with_setup
try:
    from nose.tools import assert_is_instance
except ImportError:
    from landlab.testing.tools import assert_is_instance
import os
import tempfile
import numpy as np

from landlab import ModelParameterDictionary
from landlab.core.model_parameter_dictionary import (MissingKeyError,
                                                     ParameterValueError)


def setup():
    from six import StringIO

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

    param_file = StringIO(_TEST_PARAM_DICT_FILE)
    param_dict = ModelParameterDictionary()
    param_dict.read_from_file(param_file)
    globals().update({
        'param_dict': param_dict,
        '_TEST_PARAM_DICT_FILE': _TEST_PARAM_DICT_FILE
    })


@with_setup(setup)
def test_read_file():
    all_keys = set([
        'FLOAT_VAL', 'INT_VAL', 'STRING_VAL', 'TRUE_BOOL_VAL',
        'FALSE_BOOL_VAL',
    ])
    param_list = set(param_dict.params())

    assert_equal(param_list, all_keys)


@with_setup(setup)
def test_read_file_name():
    (prm_fd, prm_file_name) = tempfile.mkstemp()

    prm_file = os.fdopen(prm_fd, 'w')
    prm_file.write(_TEST_PARAM_DICT_FILE)
    prm_file.close()

    param_dict = ModelParameterDictionary()
    param_dict.read_from_file(prm_file_name)

    os.remove(prm_file_name)

    all_keys = set([
        'FLOAT_VAL', 'INT_VAL', 'STRING_VAL', 'TRUE_BOOL_VAL',
        'FALSE_BOOL_VAL',
    ])
    param_list = set(param_dict.params())

    assert_equal(param_list, all_keys)


@with_setup(setup)
def test_read_file_like_twice():
    from six import StringIO
    param_file = StringIO(_TEST_PARAM_DICT_FILE)
    param_dict_1 = ModelParameterDictionary()
    param_dict_2 = ModelParameterDictionary()

    param_dict_1.read_from_file(param_file)
    param_dict_2.read_from_file(param_file)


@with_setup(setup)
def test_read_int():
    assert_equal(param_dict.read_int('INT_VAL'), 1)

    assert_raises(ParameterValueError, param_dict.read_int, 'FLOAT_VAL')

    assert_raises(MissingKeyError, param_dict.read_int, 'MISSING_INT')

    assert_equal(param_dict.read_int('MISSING_INT', 2), 2)


@with_setup(setup)
def test_get_int():
    assert_equal(param_dict.get('INT_VAL', ptype=int), 1)

    assert_raises(ParameterValueError, param_dict.get, 'FLOAT_VAL', ptype=int)

    assert_raises(MissingKeyError, param_dict.get, 'MISSING_INT', ptype=int)


@with_setup(setup)
def test_set_default():
    param_dict.setdefault('MISSING_INT', 2)
    assert_equal(param_dict.read_int('MISSING_INT'), 2)


@with_setup(setup)
def test_read_float():
    assert_equal(param_dict.read_float('FLOAT_VAL'), 2.2)
    assert_equal(param_dict.read_float('INT_VAL'), 1)

    assert_raises(ParameterValueError, param_dict.read_float, 'STRING_VAL')

    assert_raises(MissingKeyError, param_dict.read_float, 'MISSING_FLOAT')


@with_setup(setup)
def test_read_string():
    assert_equal(param_dict.read_string('STRING_VAL'), 'The Landlab')
    assert_equal(param_dict.read_string('INT_VAL'), '1')
    assert_equal(param_dict.read_string('FLOAT_VAL'), '2.2')

    assert_raises(MissingKeyError, param_dict.read_string, 'MISSING_STRING')


@with_setup(setup)
def test_read_bool():
    assert_equal(param_dict.read_bool('TRUE_BOOL_VAL'), True)
    assert_equal(param_dict.read_bool('FALSE_BOOL_VAL'), False)

    assert_raises(MissingKeyError, param_dict.read_bool, 'MISSING_BOOLEAN')

    assert_raises(ParameterValueError, param_dict.read_bool, 'STRING_VAL')


@with_setup(setup)
def test_dict_keys():
    all_keys = set([
        'FLOAT_VAL', 'INT_VAL', 'STRING_VAL', 'TRUE_BOOL_VAL',
        'FALSE_BOOL_VAL',
    ])
    assert_equal(set(param_dict), all_keys)
    for key in all_keys:
        assert_true(key in param_dict)


@with_setup(setup)
def test_dict_index():
    assert_equal(param_dict['INT_VAL'], '1')


def setup_auto_type():
    from six import StringIO

    _TEST_FILE = u"""
# A Comment
INT_VAL:
1
DBL_VAL:
1.2
STR_VAL:
landlab
BOOL_VAL:
true
INT_ARRAY_VAL:
1,2 ,4 ,7

DBL_ARRAY_VAL:
1.,2. ,4. ,7.
    """
    param_dict = ModelParameterDictionary(auto_type=True,
                                          from_file=StringIO(_TEST_FILE))
    globals().update({
        'param_dict': param_dict,
    })


@with_setup(setup_auto_type)
def test_auto_type():
    assert_equal(param_dict['INT_VAL'], 1)
    assert_equal(param_dict['DBL_VAL'], 1.2)
    assert_equal(param_dict['STR_VAL'], 'landlab')
    assert_equal(param_dict['BOOL_VAL'], True)


@with_setup(setup_auto_type)
def test_int_vector():
    val = param_dict['INT_ARRAY_VAL']

    assert_equal(list(val), [1, 2, 4, 7])
    assert_is_instance(val, np.ndarray)
    assert_equal(val.dtype, np.int)


@with_setup(setup_auto_type)
def test_float_vector():
    val = param_dict['DBL_ARRAY_VAL']

    assert_equal(list(val), [1., 2., 4., 7.])
    assert_is_instance(val, np.ndarray)
    assert_equal(val.dtype, np.float)

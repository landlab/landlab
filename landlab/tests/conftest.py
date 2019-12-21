import pytest
from six import StringIO

from landlab import ModelParameterDictionary


@pytest.fixture
def pdict_setup():
    test_param_dict_file = u"""
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

    param_file = StringIO(test_param_dict_file)
    param_dict = ModelParameterDictionary()
    param_dict.read_from_file(param_file)

    class ParamDictSetup(object):
        pass

    setup = ParamDictSetup()
    setup.param_dict = param_dict
    setup.param_dict_file = test_param_dict_file

    return setup


@pytest.fixture
def auto_type_setup():
    TEST_FILE = u"""
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
    param_dict = ModelParameterDictionary(auto_type=True, from_file=StringIO(TEST_FILE))

    return param_dict

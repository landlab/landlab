import pytest

from landlab import ModelParameterDictionary


@pytest.fixture
def pdict_setup():
    from six import StringIO

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

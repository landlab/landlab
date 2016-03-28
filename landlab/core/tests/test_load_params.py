from nose.tools import assert_true, assert_false, assert_raises
try:
    from nose.tools import assert_is_instance, assert_dict_equal
except ImportError:
    from landlab.testing.tools import assert_is_instance, assert_dict_equal
from six import StringIO

from landlab.core import load_params
from landlab.testing.tools import cdtemp


YAML_PARAMS_STR = """
x: 1e7
y: 1
z: [1, 2]
a: frog
"""
MPD_PARAMS_STR = """
x: A float
1e7
y: An int
1
a: A string
frog
"""
YAML_PARAMS = {
    'x': 1e7,
    'y': 1,
    'z': [1, 2],
    'a': 'frog',
}
MPD_PARAMS = {
    'x': 1e7,
    'y': 1,
    'a': 'frog',
}


def test_from_yaml_string():
    """Load parameters from YAML-formatted string."""
    params = load_params(YAML_PARAMS_STR)
    assert_dict_equal(params, YAML_PARAMS)
    assert_is_instance(params['x'], float)
    assert_is_instance(params['y'], int)


def test_from_yaml_file_like():
    """Load parameters from YAML-formatted string."""
    params = load_params(StringIO(YAML_PARAMS_STR))
    assert_dict_equal(params, YAML_PARAMS)
    assert_is_instance(params['x'], float)
    assert_is_instance(params['y'], int)


def test_from_yaml_path():
    """Load parameters from YAML-formatted string."""
    with cdtemp() as dir:
        with open('params.yaml', 'w') as fp:
            fp.write(YAML_PARAMS_STR)
        params = load_params('./params.yaml')
    assert_dict_equal(params, YAML_PARAMS)
    assert_is_instance(params['x'], float)
    assert_is_instance(params['y'], int)


def test_from_mpd_string():
    """Load parameters from YAML-formatted string."""
    params = load_params(MPD_PARAMS_STR)
    assert_dict_equal(params, MPD_PARAMS)
    assert_is_instance(params['x'], float)
    assert_is_instance(params['y'], int)


def test_from_yaml_file_like():
    """Load parameters from YAML-formatted string."""
    params = load_params(StringIO(MPD_PARAMS_STR))
    assert_dict_equal(params, MPD_PARAMS)
    assert_is_instance(params['x'], float)
    assert_is_instance(params['y'], int)


def test_from_yaml_path():
    """Load parameters from YAML-formatted string."""
    with cdtemp() as dir:
        with open('params.txt', 'w') as fp:
            fp.write(MPD_PARAMS_STR)
        params = load_params('./params.txt')
    assert_dict_equal(params, MPD_PARAMS)
    assert_is_instance(params['x'], float)
    assert_is_instance(params['y'], int)

from six import StringIO

from landlab.core import load_params

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
YAML_PARAMS = {"x": 1e7, "y": 1, "z": [1, 2], "a": "frog"}
MPD_PARAMS = {"x": 1e7, "y": 1, "a": "frog"}


def test_from_yaml_string():
    """Load parameters from YAML-formatted string."""
    params = load_params(YAML_PARAMS_STR)
    assert params == YAML_PARAMS
    assert isinstance(params["x"], float)
    assert isinstance(params["y"], int)


def test_from_yaml_file_like():
    """Load parameters from YAML-formatted string."""
    params = load_params(StringIO(YAML_PARAMS_STR))
    assert params == YAML_PARAMS
    assert isinstance(params["x"], float)
    assert isinstance(params["y"], int)


def test_from_yaml_path(tmpdir):
    """Load parameters from YAML-formatted string."""
    with tmpdir.as_cwd():
        with open("params.yaml", "w") as fp:
            fp.write(YAML_PARAMS_STR)
        params = load_params("./params.yaml")
    assert params == YAML_PARAMS
    assert isinstance(params["x"], float)
    assert isinstance(params["y"], int)


def test_from_mpd_string():
    """Load parameters from YAML-formatted string."""
    params = load_params(MPD_PARAMS_STR)
    assert params == MPD_PARAMS
    assert isinstance(params["x"], float)
    assert isinstance(params["y"], int)


def test_from_mpd_file_like():
    """Load parameters from YAML-formatted string."""
    params = load_params(StringIO(MPD_PARAMS_STR))
    assert params == MPD_PARAMS
    assert isinstance(params["x"], float)
    assert isinstance(params["y"], int)


def test_from_mpd_path(tmpdir):
    """Load parameters from YAML-formatted string."""
    with tmpdir.as_cwd():
        with open("params.txt", "w") as fp:
            fp.write(MPD_PARAMS_STR)
        params = load_params("./params.txt")
    assert params == MPD_PARAMS
    assert isinstance(params["x"], float)
    assert isinstance(params["y"], int)

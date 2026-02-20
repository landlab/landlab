from io import StringIO

import pytest

from landlab.core import load_params

YAML_PARAMS_STR = """
x: 1e7
y: 1
z: [1, 2]
a: frog
"""

NOT_A_YAML_PARAMS_STR = """
x: 1e7
y: 1
z: [1, 2]
a: frog
a
"""

YAML_PARAMS_LIST_STR = """[2,3,2,1]"""

YAML_PARAMS_EMPTY_STR = """ """

YAML_PARAMS = {"x": 1e7, "y": 1, "z": [1, 2], "a": "frog"}


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


def test_read_text_file(tmpdir):
    """load parameters from a text file"""
    with tmpdir.as_cwd():
        with open("params.txt", "w") as fp:
            fp.write(YAML_PARAMS_STR)
        params = load_params("./params.txt")
    assert params == YAML_PARAMS
    assert isinstance(params["x"], float)
    assert isinstance(params["y"], int)


def test_read_file_no_extension(tmpdir):
    """load parameters from a text file that lacks an extension"""
    with tmpdir.as_cwd():
        with open("params", "w") as fp:
            fp.write(YAML_PARAMS_STR)
        params = load_params("./params")
    assert params == YAML_PARAMS
    assert isinstance(params["x"], float)
    assert isinstance(params["y"], int)


def test_file_found_but_faulty_yaml_syntax(tmpdir):
    with tmpdir.as_cwd():
        with open("params_not_yaml.txt", "w") as fp:
            fp.write(NOT_A_YAML_PARAMS_STR)
        with pytest.raises(ValueError):
            load_params("./params_not_yaml.txt")


def test_file_found_but_not_a_dict(tmpdir):
    with tmpdir.as_cwd():
        with open("params_list.txt", "w") as fp:
            fp.write(YAML_PARAMS_LIST_STR)
        with pytest.raises(
            ValueError,
            match="File found but parsing produced a list instead of a dictionary. Ensure that your file uses 'key: value' syntax",
        ):
            load_params("./params_list.txt")


def test_yaml_string_but_is_not_a_dict():
    yaml_string = YAML_PARAMS_LIST_STR
    with pytest.raises(
        ValueError,
        match="The YAML content was parsed successfully but produced a list instead of a dictionary. Ensure that your file uses 'key: value' syntax",
    ):
        load_params(yaml_string)


def test_yaml_string_is_empty():
    yaml_string = """ """
    with pytest.raises(ValueError, match="The parameter file is empty"):
        load_params(yaml_string)


def test_file_found_but_is_empty(tmpdir):
    with tmpdir.as_cwd():
        with open("params_empty.txt", "w") as fp:
            fp.write(YAML_PARAMS_EMPTY_STR)
        with pytest.raises(ValueError, match="File found but is empty"):
            load_params("./params_empty.txt")

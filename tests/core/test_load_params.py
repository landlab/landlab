from contextlib import chdir
from io import BytesIO
from io import StringIO
from unittest.mock import patch

import pytest
import yaml

from landlab.core.model_parameter_loader import EmptyConfigurationError
from landlab.core.model_parameter_loader import InvalidConfigurationError
from landlab.core.model_parameter_loader import load_file_contents
from landlab.core.model_parameter_loader import load_params

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
        with pytest.raises(ValueError, match="^params_list.txt:"):
            load_params("params_list.txt")


def test_yaml_string_but_is_not_a_dict():
    yaml_string = YAML_PARAMS_LIST_STR
    with pytest.raises(ValueError, match="^expected a parameter dictionary"):
        load_params(yaml_string)


def test_yaml_string_is_empty():
    yaml_string = """ """
    with pytest.raises(ValueError, match="^expected a parameter dictionary"):
        load_params(yaml_string)


def test_file_found_but_is_empty(tmpdir):
    with tmpdir.as_cwd():
        with open("params_empty.txt", "w") as fp:
            fp.write(YAML_PARAMS_EMPTY_STR)
        with pytest.raises(
            ValueError, match="^params_empty.txt: expected a parameter dictionary"
        ):
            load_params("params_empty.txt")


@pytest.mark.parametrize(
    "contents, expected",
    [
        ("x: 1.5", {"x": 1.5}),
        ("start: 0.\nstop: 10.\nstep: 2.\n", {"start": 0.0, "stop": 10.0, "step": 2.0}),
        ("output: results.nc", {"output": "results.nc"}),
        ("x: " + "a" * 4096, {"x": "a" * 4096}),
        ("{}", {}),
    ],
    ids=["decimal", "multiline-decimals", "filename-value", "long-text", "empty-dict"],
)
def test_inline_yaml_edge_cases(contents, expected):
    assert load_params(contents) == expected


@pytest.mark.parametrize(
    "contents", ("{x: 1.5}", "start: 0.0\nstop: 10.0", "output: foobar.txt")
)
def test_inline_yaml_decimals_are_not_extensions(contents):
    """Check decimals not interpreted as file extensions."""
    assert load_params(contents) == yaml.safe_load(contents)


def test_from_pathlib_path(tmp_path):
    path = tmp_path / "params.yaml"
    path.write_text("foo: bar")
    assert load_params(path) == {"foo": "bar"}


def test_from_binary_stream():
    assert load_params(BytesIO(b"foo: bar")) == {"foo": "bar"}


@pytest.mark.parametrize("contents", ("", " \n", "# only a comment\n", "null", "~"))
def test_empty_configuration(contents):
    with pytest.raises(EmptyConfigurationError, match="YAML content was empty$"):
        load_params(contents)


@pytest.mark.parametrize("contents", ["[1, 2]", "123", "true", "a scalar"])
def test_non_dictionary_stream(contents):
    with pytest.raises(
        InvalidConfigurationError, match="expected a parameter dictionary"
    ):
        load_params(StringIO(contents))


def test_malformed_yaml():
    with pytest.raises(
        InvalidConfigurationError, match="could not parse parameters as YAML"
    ):
        load_params("x: [")


@pytest.mark.parametrize("filename", ["missing.yaml", "missing", "123"])
def test_missing_string_filename(filename, tmp_path):
    with chdir(tmp_path):
        with pytest.raises(
            InvalidConfigurationError, match="check that the file exists"
        ):
            load_params(filename)


def test_missing_path(tmp_path):
    with pytest.raises(FileNotFoundError):
        load_params(tmp_path / "missing.yaml")


@pytest.mark.parametrize("value", [None, 42, [], {}])
def test_bad_input_type(value):
    with pytest.raises(TypeError, match="^'file_like' must be"):
        load_params(value)


def test_keys_must_be_strings():
    with pytest.raises(
        InvalidConfigurationError, match="^parameter names must be strings."
    ):
        load_params("1: foo")


def test_oserror_propagates():
    with patch("builtins.open") as mock_open:
        mock_open.side_effect = PermissionError("secret.yaml")
        with pytest.raises(PermissionError, match="^secret.yaml"):
            load_params("secret.yaml")


@pytest.mark.parametrize("contents", ("foobar", "foo\nbar"))
def test_load_file_contents_matches(contents):
    with patch("landlab.core.model_parameter_loader._read_parameter_source") as reader:
        reader.return_value = (contents, "text")
        assert load_file_contents(contents) == contents
    reader.assert_called_once()

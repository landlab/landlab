import errno
import os
import pathlib
import re
import tomllib
from collections.abc import Callable
from typing import Any
from typing import BinaryIO
from typing import TextIO

import yaml
from requireit import require_one_of

_loader = yaml.SafeLoader
_loader.add_implicit_resolver(
    "tag:yaml.org,2002:float",
    re.compile(
        """^(?:
               [-+]?(?:[0-9][0-9_]*)\\.[0-9_]*(?:[eE][-+]?[0-9]+)?
               |[-+]?(?:[0-9][0-9_]*)(?:[eE][-+]?[0-9]+)
               |\\.[0-9_]+(?:[eE][-+][0-9]+)?
               |[-+]?[0-9][0-9_]*(?::[0-5]?[0-9])+\\.[0-9_]*
               |[-+]?\\.(?:inf|Inf|INF)
               |\\.(?:nan|NaN|NAN))$""",
        re.X,
    ),
    list("-+0123456789."),
)


type ParameterLoader = Callable[[str], dict[str, Any]]


class ConfigurationError(ValueError):
    pass


class InvalidConfigurationError(ConfigurationError):
    pass


class NonDictConfigurationError(InvalidConfigurationError):
    pass


def load_file_contents(
    file_like: str | pathlib.Path | TextIO | BinaryIO,
) -> str:
    """Load the contents of a file or file-like object.

    Parameters
    ----------
    file_like : file_like or str
        File to load either as a file-like object, path to an existing file,
        or the contents of a file.

    Returns
    -------
    str
        The contents of the file.
    """
    contents, _ = _read_parameter_source(file_like)
    return contents


def _read_parameter_source(
    file_like: str | pathlib.Path | TextIO | BinaryIO,
) -> tuple[str, str]:
    if hasattr(file_like, "read"):
        contents = file_like.read()
        if isinstance(contents, bytes):
            contents = contents.decode("utf-8")
        return contents, "stream"

    if isinstance(file_like, pathlib.Path):
        with open(file_like) as fp:
            return fp.read(), "file"

    if isinstance(file_like, str):
        try:
            with open(file_like) as fp:
                return fp.read(), "file"
        except OSError as error:
            if error.errno in (
                errno.EINVAL,
                errno.ENAMETOOLONG,
                errno.ENOENT,
            ):
                return file_like, "text"
            raise

    raise TypeError(
        "'file_like' must be either a file-like object, a path to a file or a string"
    )


def _load_yaml(contents: str) -> dict[str, Any]:
    try:
        params = yaml.load(contents, Loader=_loader)
    except yaml.YAMLError as exc:
        raise InvalidConfigurationError(
            f"could not parse parameters as YAML: {exc}"
        ) from exc

    if params is None:
        params = {}

    if not isinstance(params, dict):
        raise NonDictConfigurationError(
            f"expected a parameter dictionary, but YAML parsing produced"
            f" {type(params).__name__}. Use 'key: value' entries."
        )

    invalid_keys = [key for key in params if not isinstance(key, str)]
    if invalid_keys:
        raise InvalidConfigurationError("parameter names must be strings.")

    return params


def _load_toml(contents: str) -> dict[str, Any]:
    try:
        return tomllib.loads(contents)
    except tomllib.TOMLDecodeError as exc:
        raise InvalidConfigurationError(
            f"could not parse parameters as TOML: {exc}"
        ) from exc


LOADERS: dict[str, ParameterLoader] = {
    ".yaml": _load_yaml,
    ".yml": _load_yaml,
    ".toml": _load_toml,
    "yaml": _load_yaml,
    "yml": _load_yaml,
    "toml": _load_toml,
}


def _get_loader(
    file_like: str | pathlib.Path | TextIO | BinaryIO,
    *,
    source: str = "text",
    fmt: str | None = None,
) -> ParameterLoader:
    if fmt is not None:
        return LOADERS[fmt]

    if source == "text":
        return LOADERS["yaml"]

    path = None
    if hasattr(file_like, "name"):
        if isinstance(file_like.name, str):
            path = file_like.name
    else:
        path = str(file_like)

    if path is not None:
        _, ext = os.path.splitext(path)
        return LOADERS.get(ext, _load_yaml)

    return LOADERS["yaml"]


def load_params(
    file_like: str | pathlib.Path | TextIO | BinaryIO,
    *,
    fmt: str | None = None,
) -> dict[str, Any]:
    """Load parameters from a YAML or TOML file.

    Parameters
    ----------
    file_like : file_like or str
        Contents of a parameter file, a file-like object, or the path to
        a parameter file.
    fmt : {'yaml', 'toml'}, optional
        Format of the parameter file. If not provided, the format is
        chosen based on the filename. Filenames, with a `.toml` extension
        are parsed as TOML, otherwise they are parsed as YAML.

    Returns
    -------
    dict
        Parameters as key-value pairs.

    Examples
    --------
    >>> from landlab.core import load_params
    >>> contents = '''
    ... start: 0.
    ... stop: 10.
    ... step: 2.
    ... '''
    >>> params = load_params(contents)
    >>> isinstance(params, dict)
    True
    >>> params["start"], params["stop"], params["step"]
    (0.0, 10.0, 2.0)
    """
    if fmt is not None:
        fmt = require_one_of(fmt, allowed=("yaml", "toml"), name="fmt")

    contents, source = _read_parameter_source(file_like)

    loader = _get_loader(file_like, source=source, fmt=fmt)
    try:
        params = loader(contents)
    except InvalidConfigurationError as error:
        msg = str(error)
        if source == "text" and isinstance(error, NonDictConfigurationError):
            msg = f"{msg} If you intended a filename, check that the file exists."
        if source == "file":
            msg = f"{file_like}: {msg}"
        raise InvalidConfigurationError(msg) from error

    return params

import errno
import pathlib
import re
from typing import Any
from typing import BinaryIO
from typing import TextIO

import yaml

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


class EmptyConfigurationError(ValueError):
    pass


class InvalidConfigurationError(ValueError):
    pass


def load_file_contents(
    file_like: str | pathlib.Path | TextIO | BinaryIO,
) -> dict[str, Any]:
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
        return file_like.read(), "stream"

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


def load_params(
    file_like: str | pathlib.Path | TextIO | BinaryIO,
) -> dict[str, Any]:
    """Load parameters from a YAML style file.

    Parameters
    ----------
    file_like : file_like or str
        Contents of a parameter file, a file-like object, or the path to
        a parameter file.

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
    contents, source = _read_parameter_source(file_like)
    prefix = f"{file_like}: " if source == "file" else ""

    try:
        params = yaml.load(contents, Loader=_loader)
    except yaml.YAMLError as exc:
        raise InvalidConfigurationError(
            f"{prefix}could not parse parameters as YAML: {exc}"
        ) from exc

    if params is None:
        raise EmptyConfigurationError(
            f"{prefix}expected a parameter dictionary, but YAML content was empty"
        )

    if not isinstance(params, dict):
        msg = (
            f"expected a parameter dictionary, but YAML parsing produced"
            f" {type(params).__name__}. Use 'key: value' entries."
        )

        if source == "text":
            msg += (
                " Input could not be opened as a file and did not contain"
                " a YAML parameter dictionary. If you intended a filename,"
                " check that the file exists."
            )

        raise InvalidConfigurationError(prefix + msg)

    invalid_keys = [key for key in params if not isinstance(key, str)]
    if invalid_keys:
        raise InvalidConfigurationError(f"{prefix}parameter names must be strings.")

    return params

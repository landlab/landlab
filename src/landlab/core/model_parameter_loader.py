import os
import re

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


def load_file_contents(file_like):
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
    try:
        contents = file_like.read()
    except AttributeError:  # was a str
        if os.path.isfile(file_like):
            with open(file_like) as fp:
                contents = fp.read()
        else:
            contents = file_like

    return contents


def load_params(file_like):
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
    contents = load_file_contents(file_like)

    params = yaml.load(contents, Loader=_loader)

    if not isinstance(params, dict):
        raise ValueError("parsing of parameter file did not produce a dict-like object")

    return params

import re

import six
import yaml

from .model_parameter_dictionary import ModelParameterDictionary


_loader = yaml.SafeLoader
_loader.add_implicit_resolver(
    u'tag:yaml.org,2002:float',
    re.compile(u'''^(?:
               [-+]?(?:[0-9][0-9_]*)\\.[0-9_]*(?:[eE][-+]?[0-9]+)?
               |[-+]?(?:[0-9][0-9_]*)(?:[eE][-+]?[0-9]+)
               |\\.[0-9_]+(?:[eE][-+][0-9]+)?
               |[-+]?[0-9][0-9_]*(?::[0-5]?[0-9])+\\.[0-9_]*
               |[-+]?\\.(?:inf|Inf|INF)
               |\\.(?:nan|NaN|NAN))$''', re.X),
    list(u'-+0123456789.'))


def load_params(file_like):
    """Load parameters from a file.

    Parameters
    ----------
    file_like : file_like or str
        Contents of a parameter file or a file-like object.

    Returns
    -------
    dict
        Parameters as key-value pairs.

    Examples
    --------
    >>> from landlab.core import load_params
    >>> contents = \"\"\"
    ... start: 0.
    ... stop: 10.
    ... step: 2.
    ... \"\"\"
    >>> params = load_params(contents)
    >>> isinstance(params, dict)
    True
    >>> params['start'], params['stop'], params['step']
    (0.0, 10.0, 2.0)

    >>> contents = \"\"\"
    ... start: Start time
    ... 0.
    ... stop: Stop time
    ... 10.
    ... step: Step time
    ... 2.
    ... \"\"\"
    >>> params = load_params(contents)
    >>> isinstance(params, dict)
    True
    >>> params['start'], params['stop'], params['step']
    (0.0, 10.0, 2.0)
    """
    try:
        params = yaml.load(file_like, Loader=_loader)
    except yaml.YAMLError:
        if isinstance(file_like, six.string_types):
            file_like = six.StringIO(file_like)
        params = ModelParameterDictionary(from_file=file_like, auto_type=True)

    if not isinstance(params, dict):
        raise ValueError(
            'parsing of parameter file did not produce a dict-like object')

    return params

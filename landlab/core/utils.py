#! /usr/bin/env python

import numpy as np


SIZEOF_INT = np.dtype(np.int).itemsize


def as_id_array(x):
    """Convert an array to an array of ids.

    Parameters
    ----------
    x : ndarray
        Array of IDs.

    Returns
    -------
    ndarray
        A, possibly new, array of IDs.

    Examples
    --------
    >>> import numpy as np
    >>> from landlab.core.utils import as_id_array
    >>> x = np.arange(5)
    >>> y = as_id_array(x)
    >>> y
    array([0, 1, 2, 3, 4])
    >>> x is y
    True

    >>> x = np.arange(5, dtype=np.int32)
    >>> y = as_id_array(x)
    >>> y
    array([0, 1, 2, 3, 4])
    >>> y.dtype == np.int
    True

    >>> x = np.arange(5, dtype=np.int64)
    >>> y = as_id_array(x)
    >>> y
    array([0, 1, 2, 3, 4])
    >>> y.dtype == np.int
    True

    >>> x = np.arange(5, dtype=np.intp)
    >>> y = as_id_array(x)
    >>> y
    array([0, 1, 2, 3, 4])
    >>> y.dtype == np.int
    True
    """
    if x.dtype == np.intp and np.dtype(np.intp) != np.int:
        return np.array(x, dtype=np.int)

    if x.dtype != np.int or x.itemsize != SIZEOF_INT:
        id_array = x.astype(np.int, copy=True)
        if id_array.dtype.itemsize != SIZEOF_INT:
            raise RuntimeError('id array is of type {dtype}'.format(dtype=id_array.dtype))
    else:
        id_array = x
    return id_array


def get_functions_from_module(mod, pattern=None):
    import inspect, re

    funcs = {}
    for name, func in inspect.getmembers(mod, inspect.isroutine):
        if pattern is None or re.match(pattern, name):
            funcs[name] = func
    return funcs


def add_functions_to_class(cls, funcs):
    for name, func in funcs.items():
        setattr(cls, name, func)


def add_module_functions_to_class(cls, module, pattern=None):
    import inspect, imp, os

    caller = inspect.stack()[1]
    path = os.path.join(os.path.dirname(caller[1]), os.path.dirname(module))

    (module, _) = os.path.splitext(os.path.basename(module))

    mod = imp.load_module(module, *imp.find_module(module, [path]))

    funcs = get_functions_from_module(mod, pattern=pattern)
    add_functions_to_class(cls, funcs)



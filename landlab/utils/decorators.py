"""General decorators for the landlab library."""

import warnings
import types
from functools import wraps

import numpy as np


class use_field_name_or_array(object):

    """Decorate a function so that it accepts a field name or array.

    Parameters
    ----------
    func : function
        A function that accepts a grid and array as arguments.

    Returns
    -------
    func
        A wrapped function that accepts a grid and either a field name or
        a numpy array.
    """

    def __init__(self, at_element):
        """Initialize the decorator.

        Parameters
        ----------
        at_element : str
            The element type that the field is defined on ('node', 'cell',
            etc.)
        """
        self._at = at_element

    def __call__(self, func):
        @wraps(func)
        def _wrapped(grid, vals, *args, **kwds):
            if isinstance(vals, types.StringTypes):
                vals = grid[self._at][vals]
            else:
                vals = np.asarray(vals).flatten()

            return func(grid, vals, *args, **kwds)
        return _wrapped


def make_return_array_immutable(func):
    """Decorate a function so that its return array is read-only.

    Parameters
    ----------
    func : function
        A function that returns a numpy array.

    Returns
    -------
    func
        A wrapped function that returns a read-only view of an array.
    """
    @wraps(func)
    def _wrapped(self, *args, **kwds):
        array = func(self, *args, **kwds)
        immutable_array = array.view()
        immutable_array.flags.writeable = False
        return immutable_array
    return _wrapped


def deprecated(func):
    """Mark a function as deprecated.

    Parameters
    ----------
    func : function
        A function.

    Returns
    -------
    func
        A wrapped function that issues a deprecation warning.
    """
    @wraps(func)
    def _wrapped(*args, **kwargs):
        warnings.warn("Call to deprecated function {}.".format(func.__name__),
                      category=DeprecationWarning)
        return func(*args, **kwargs)
    _wrapped.__dict__.update(func.__dict__)

    return _wrapped

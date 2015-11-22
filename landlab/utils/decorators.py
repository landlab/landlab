"""General decorators for the landlab library."""

import warnings
from functools import wraps

import numpy as np
import six


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

    Examples
    --------
    >>> from landlab import RasterModelGrid
    >>> grid = RasterModelGrid((4, 5), spacing=(1, 2))

    >>> def my_func(grid, vals):
    ...     return grid.area_of_cell * vals
    >>> my_func(grid, np.arange(grid.number_of_cells))
    array([  0.,   2.,   4.,   6.,   8.,  10.])

    Decorate the function so that the second argument can be array-like or
    the name of a field contained withing the grid. The decorator takes a
    single argument that is the name (as a `str`) of the grid element that
    the values are defined on ("node", "cell", etc.).

    >>> from landlab.utils.decorators import use_field_name_or_array
    >>> @use_field_name_or_array('cell')
    ... def my_func(grid, vals):
    ...     return grid.area_of_cell * vals

    The array of values now can be list or anything that can be converted to
    a numpy array.

    >>> my_func(grid, [0, 1, 2, 3, 4, 5])
    array([  0.,   2.,   4.,   6.,   8.,  10.])

    The array of values doesn't have to be flat.

    >>> vals = np.array([[0, 1, 2], [3, 4, 5]])
    >>> my_func(grid, vals)
    array([  0.,   2.,   4.,   6.,   8.,  10.])

    The array of values can be a field name.

    >>> _ = grid.add_field('cell', 'elevation', [0, 1, 2, 3, 4, 5])
    >>> my_func(grid, 'elevation')
    array([  0.,   2.,   4.,   6.,   8.,  10.])
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
        """Wrap the function."""
        @wraps(func)
        def _wrapped(grid, vals, *args, **kwds):
            """Convert the second argument to an array."""
            if isinstance(vals, six.string_types):
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
        """Make the returned array read-only."""
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
        """Warn that the function is deprecated before calling it."""
        warnings.warn(
            "Call to deprecated function {name}.".format(name=func.__name__),
            category=DeprecationWarning)
        return func(*args, **kwargs)
    _wrapped.__dict__.update(func.__dict__)

    return _wrapped

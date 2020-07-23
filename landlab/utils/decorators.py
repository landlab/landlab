"""General decorators for the landlab library.

General Landlab decorators
++++++++++++++++++++++++++

.. autosummary::

    ~landlab.utils.decorators.use_field_name_or_array
    ~landlab.utils.decorators.make_return_array_immutable
    ~landlab.utils.decorators.deprecated
"""

import inspect
import os
import textwrap
import warnings
from functools import wraps

import numpy as np

from landlab import FieldError


class cache_result_in_object(object):
    def __init__(self, cache_as=None):
        self._attr = cache_as

    def __call__(self, func):
        name = self._attr or "_" + func.__name__

        @wraps(func)
        def _wrapped(obj):
            if not hasattr(obj, name):
                setattr(obj, name, func(obj))
            return getattr(obj, name)

        return _wrapped


class store_result_in_grid(object):
    def __init__(self, name=None):
        self._attr = name

    def __call__(self, func):
        @wraps(func)
        def _wrapped(grid):
            name = self._attr or "_" + func.__name__
            try:
                getattr(grid, name)
            except AttributeError:
                setattr(grid, name, func(grid))
            finally:
                return getattr(grid, name)

        return _wrapped


class store_result_in_dataset(object):
    def __init__(self, dataset=None, name=None):
        self._dataset = dataset
        self._attr = name

    def __call__(self, func):
        @wraps(func)
        def _wrapped(grid):
            name = self._attr or func.__name__
            if self._dataset:
                ds = getattr(grid, self._dataset)
            else:
                ds = grid

            if name not in ds:
                setattr(grid, self._dataset, ds.update(dict(name=func(grid))))

            return getattr(grid, self._dataset)[name].values

        return _wrapped


def read_only_array(func):
    """Decorate a function so that its return array is read-only.

    Parameters
    ----------
    func : function
        A function that returns a numpy array.

    Returns
    -------
    func
        A wrapped function that returns a read-only numpy array.
    """

    @wraps(func)
    def _wrapped(self, *args, **kwds):
        """Make the returned array read-only."""
        array = func(self, *args, **kwds)
        array.flags.writeable = False
        return array

    return _wrapped


def add_signature_to_doc(func):
    """Add a function signature to the first line of a docstring.

    Add the signature of a function to the first line of its docstring.
    This is useful for functions that are decorated in such a way
    that its signature changes.

    Parameters
    ----------
    func : function
        Function to add signature to.

    Returns
    -------
    str
        The new docstring with the function signature on the first line.

    Examples
    --------
    >>> from landlab.utils.decorators import add_signature_to_doc

    >>> def foo(arg1, kwd=None):
    ...     '''Do something.'''
    ...     pass
    >>> print(add_signature_to_doc(foo))
    foo(arg1, kwd=None)
    <BLANKLINE>
    Do something.
    """
    return """{name}{argspec}

{body}""".format(
        name=func.__name__,
        argspec=str(inspect.signature(func)),
        body=inspect.getdoc(func),
    )


class use_field_name_or_array(object):

    """Decorate a function so that it accepts a field name or array.

    Parameters
    ----------
    func : function
        A function that accepts a grid and array as arguments.
    at_element : str
        The element type that the field is defined on ('node', 'cell',
        etc.)

    Returns
    -------
    func
        A wrapped function that accepts a grid and either a field name or
        a numpy array.

    Examples
    --------
    >>> from landlab import RasterModelGrid
    >>> grid = RasterModelGrid((4, 5), xy_spacing=(1, 2))

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

    >>> _ = grid.add_field("elevation", [0, 1, 2, 3, 4, 5], at="cell")
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

        func.__doc__ = add_signature_to_doc(func)

        @wraps(func)
        def _wrapped(grid, vals, *args, **kwds):
            """Convert the second argument to an array."""
            if isinstance(vals, str):
                if vals in grid[self._at]:
                    vals = grid[self._at][vals]
                else:
                    raise FieldError(vals)
            else:
                vals = np.asarray(vals).flatten()

            return func(grid, vals, *args, **kwds)

        return _wrapped


class use_field_name_array_or_value(object):

    """Decorate a function so that it accepts a field name, array, or value.

    Parameters
    ----------
    func : function
        A function that accepts a grid and array as arguments.
    at_element : str
        The element type that the field is defined on ('node', 'cell',
        etc.)

    Returns
    -------
    func
        A wrapped function that accepts a grid and either a field name,
        a numpy array, or a value as arguments.

    Examples
    --------
    >>> from landlab import RasterModelGrid
    >>> grid = RasterModelGrid((4, 5), xy_spacing=(1, 2))

    >>> def my_func(grid, vals):
    ...     return grid.area_of_cell * vals
    >>> my_func(grid, np.arange(grid.number_of_cells))
    array([  0.,   2.,   4.,   6.,   8.,  10.])

    Decorate the function so that the second argument can be array-like,
    the name of a field contained withing the grid, or a value (float, int,
    etc.). The decorator takes a single argument that is the name (as a `str`)
    of the grid element that the values are defined on ("node", "cell", etc.).

    >>> from landlab.utils.decorators import use_field_name_array_or_value
    >>> @use_field_name_array_or_value('cell')
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

    >>> _ = grid.add_field("elevation", [0, 1, 2, 3, 4, 5], at="cell")
    >>> my_func(grid, 'elevation')
    array([  0.,   2.,   4.,   6.,   8.,  10.])

    The array of values can be a value (float, int, etc.).

    >>> my_func(grid, 4.0)
    array([ 8.,  8.,  8.,  8.,  8.,  8.])
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

        func.__doc__ = add_signature_to_doc(func)

        @wraps(func)
        def _wrapped(grid, vals, *args, **kwds):
            """Convert the second argument to an array."""
            if isinstance(vals, str):
                if vals in grid[self._at]:
                    vals = grid[self._at][vals]
                else:
                    raise FieldError(vals)
            else:
                expected_size = grid.size(self._at)
                vals = np.asarray(vals).ravel()
                if vals.size == 1:
                    vals = np.broadcast_to(vals, (expected_size,))

                if vals.size != expected_size:
                    raise ValueError(
                        (
                            "Array passed to function decorated with "
                            "use_field_name_array_or_value is not "
                            "the size of fields at " + self._at
                        )
                    )
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


def deprecated(use, version):
    """Mark a function as deprecated.

    Parameters
    ----------
    use : str
        Name of replacement function to use.
    version : str
        Version number when function was marked as deprecated.

    Returns
    -------
    func
        A wrapped function that issues a deprecation warning.
    """

    def real_decorator(func):
        warning_str = """
.. note:: This method is deprecated as of Landlab version {ver}.

    Use :func:`{use}` instead.

""".format(
            ver=version, use=use
        )

        doc_lines = (func.__doc__ or "").split(os.linesep)

        for lineno, line in enumerate(doc_lines):
            if len(line.rstrip()) == 0:
                break

        head = doc_lines[:lineno]
        body = doc_lines[lineno:]

        head = textwrap.dedent(os.linesep.join(head))
        body = textwrap.dedent(os.linesep.join(body))

        func.__doc__ = os.linesep.join([head, warning_str, body])

        @wraps(func)
        def _wrapped(*args, **kwargs):
            if func.__name__.startswith("_"):
                pass
            else:
                warnings.warn(
                    message="Call to deprecated function {name}.".format(
                        name=func.__name__
                    ),
                    category=DeprecationWarning,
                )
            return func(*args, **kwargs)

        _wrapped.__dict__.update(func.__dict__)

        return _wrapped

    return real_decorator

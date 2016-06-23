"""General decorators for the landlab library.

General Landlab decorators
++++++++++++++++++++++++++

.. autosummary::
    :toctree: generated/

    ~landlab.utils.decorators.use_file_name_or_kwds
    ~landlab.utils.decorators.use_field_name_or_array
    ~landlab.utils.decorators.make_return_array_immutable
    ~landlab.utils.decorators.deprecated
"""

import os
import warnings
from functools import wraps
import textwrap

import numpy as np
import six

from ..core.model_parameter_loader import load_params


def use_file_name_or_kwds(func):

    """Decorate a method so that it takes a file name or keywords.

    Parameters
    ----------
    func : A function
        A method function that accepts a ModelGrid as a first argument.

    Returns
    -------
    function
        A function that takes an optional second argument, a file, from which
        to read keywords.

    Examples
    --------
    >>> from landlab import RasterModelGrid
    >>> from landlab.utils.decorators import use_file_name_or_kwds

    >>> class MyClass(object):
    ...     @use_file_name_or_kwds
    ...     def __init__(self, grid, kw=0.):
    ...         self.kw = kw

    >>> grid = RasterModelGrid((4, 5))
    >>> foo = MyClass(grid)
    >>> foo.kw
    0.0

    >>> foo = MyClass(grid, "kw: 1945")
    >>> foo.kw
    1945
    >>> foo = MyClass(grid, "kw: 1945", kw=1973)
    >>> foo.kw
    1973

    >>> mpd = \"\"\"
    ... kw: kw value
    ... 1e6
    ... \"\"\"
    >>> foo = MyClass(grid, mpd)
    >>> foo.kw
    1000000.0
    """

    @wraps(func)
    def _wrapped(self, *args, **kwds):
        from ..grid import ModelGrid

        if not isinstance(args[0], ModelGrid):
            raise ValueError('first argument must be a ModelGrid')

        if len(args) == 2:
            warnings.warn(
                "Passing a file to a component's __init__ method is "
                "deprecated. Instead, pass parameters as keywords.",
                category=DeprecationWarning)

            if os.path.isfile(args[1]):
                with open(args[1], 'r') as fp:
                    params = load_params(fp)
            else:
                params = load_params(args[1])
        else:
            params = {}

        params.update(kwds)

        func(self, args[0], **params)

    return _wrapped


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


# def deprecated(use, version):
#     """Mark a function as deprecated.
# 
#     Parameters
#     ----------
#     use : str
#         Name of replacement function to use.
#     version : str
#         Version number when function was marked as deprecated.
# 
#     Returns
#     -------
#     func
#         A wrapped function that issues a deprecation warning.
#     """
#     def real_decorator(func):
# 
#         def _wrapped(*args, **kwargs):
#             """Warn that the function is deprecated before calling it."""
#             warnings.warn(
#                 "Call to deprecated function {name}.".format(
#                     name=func.__name__), category=DeprecationWarning)
#             return func(*args, **kwargs)
#         _wrapped.__dict__.update(func.__dict__)
# 
#         return _wrapped
# 
#     return real_decorator

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

""".format(ver=version, use=use)

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
            if func.__name__.startswith('_'):
                pass
            else:
                warnings.warn(
                    message="Call to deprecated function {name}.".format(
                        name=func.__name__), category=DeprecationWarning)
            return func(*args, **kwargs)
        _wrapped.__dict__.update(func.__dict__)

        return _wrapped

    return real_decorator

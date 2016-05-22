#! /usr/bin/env python
"""Some utilities for the landlab package.

Landlab utilities
+++++++++++++++++

.. autosummary::
    :toctree: generated/

    ~landlab.core.utils.radians_to_degrees
    ~landlab.core.utils.extend_array
    ~landlab.core.utils.as_id_array
    ~landlab.core.utils.make_optional_arg_into_id_array
    ~landlab.core.utils.get_functions_from_module
    ~landlab.core.utils.add_functions_to_class
    ~landlab.core.utils.add_module_functions_to_class
    ~landlab.core.utils.strip_grid_from_method_docstring
    ~landlab.core.utils.argsort_points_by_x_then_y
    ~landlab.core.utils.sort_points_by_x_then_y
    ~landlab.core.utils.anticlockwise_argsort_points
"""

import numpy as np



SIZEOF_INT = np.dtype(np.int).itemsize


def radians_to_degrees(rads):
    """Convert radians to compass-style degrees.

    Convert angles (measured counter-clockwise from the positive x-axis) in
    radians to angles in degrees measured clockwise starting from north.

    Parameters
    ----------
    rads : float or ndarray
        Angles in radians.

    Returns
    -------
    degrees : float or ndarray
        Converted angles in degrees.

    Examples
    --------
    >>> import numpy as np
    >>> from landlab.core.utils import radians_to_degrees

    >>> radians_to_degrees(0.)
    90.0
    >>> radians_to_degrees(np.pi / 2.)
    0.0
    >>> radians_to_degrees(- 3 * np.pi / 2.)
    0.0
    >>> radians_to_degrees(np.array([- np.pi, np.pi]))
    array([ 270.,  270.])
    """
    degrees = (5. * np.pi / 2. - rads) % (2. * np.pi)
    return 180. / np.pi * degrees


def extend_array(x, fill=0):
    """Extend an array by one element.

    The new array will appear as a view of the input array. However, its
    data now points to a newly-allocated buffer that's one element longer
    and contains a copy of the contents of *x*. The last element of the
    buffer is filled with *fill*. To access the extended array, use the
    *x* attribute of the returned array.

    Parameters
    ----------
    x : ndarray
        The array to extend.
    fill : number, optional
        The value to fill the extra element with.

    Returns
    -------
    ndarray
        A view of the original array with the data buffer extended by one
        element.

    Examples
    --------
    >>> from landlab.core.utils import extend_array
    >>> import numpy as np
    >>> arr = extend_array(np.arange(4).reshape((2, 2)))
    >>> arr
    array([[0, 1],
           [2, 3]])
    >>> arr.ext
    array([0, 1, 2, 3, 0])

    If the array is already extended, don't extend it more. However,
    possibly change its fill value.

    >>> rtn = extend_array(arr, fill=999)
    >>> rtn is arr
    True
    >>> rtn.ext
    array([  0,   1,   2,   3, 999])
    """
    if hasattr(x, 'ext'):
        x.ext[-1] = fill
        return x

    extended = np.empty(x.size + 1, dtype=x.dtype)
    extended[:-1] = x.flat
    extended[-1] = fill
    x.data = extended.data

    class array(np.ndarray):
        def __new__(cls, arr):
            """Instantiate the class with a view of the base array."""
            obj = np.asarray(arr).view(cls)
            obj.ext = extended
            return obj

    return array(x)


def as_id_array(array):
    """Convert an array to an array of ids.

    Parameters
    ----------
    array : ndarray
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

    >>> x = np.arange(5, dtype=np.int)
    >>> y = as_id_array(x)
    >>> y
    array([0, 1, 2, 3, 4])

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
    if array.dtype == np.int:
        return array.view(np.int)
    else:
        return array.astype(np.int)


if np.dtype(np.intp) == np.int:
    def _as_id_array(array):
        if array.dtype == np.intp or array.dtype == np.int:
            return array.view(np.int)
        else:
            return array.astype(np.int)
else:
    def _as_id_array(array):
        if array.dtype == np.int:
            return array.view(np.int)
        else:
            return array.astype(np.int)


def make_optional_arg_into_id_array(number_of_elements, *args):
    """Transform an optional argument into an array of element ids.

    Many landlab functions an optional argument of element ids that tells the
    function to operate only on the elements provided. However, if the argument
    is absent, all of the elements are to be operated on. This is a convenience
    function that converts such an arguments list into an array of elements
    ids.

    Parameters
    ----------
    number_of_elements : int
        Number of elements in the grid.
    array : array_like
        Iterable to convert to an array.

    Returns
    -------
    ndarray
        Input array converted to a numpy array, or a newly-created numpy
        array.

    Examples
    --------
    >>> import numpy as np
    >>> from landlab.core.utils import make_optional_arg_into_id_array
    >>> make_optional_arg_into_id_array(4)
    array([0, 1, 2, 3])
    >>> make_optional_arg_into_id_array(4, [0, 0, 0, 0])
    array([0, 0, 0, 0])
    >>> make_optional_arg_into_id_array(4, (1, 1, 1, 1))
    array([1, 1, 1, 1])
    >>> make_optional_arg_into_id_array(4, np.ones(4))
    array([1, 1, 1, 1])
    >>> make_optional_arg_into_id_array(4, 0)
    array([0])
    >>> make_optional_arg_into_id_array(4, np.array([[1, 2], [3, 4]]))
    array([1, 2, 3, 4])
    """
    if len(args) == 0:
        ids = np.arange(number_of_elements, dtype=np.int)
    elif len(args) == 1:
        ids = as_id_array(np.asarray(args[0])).reshape((-1, ))
    else:
        raise ValueError('Number of arguments must be 0 or 1.')

    return ids


def get_functions_from_module(mod, pattern=None):
    """Get all the function in a module.

    Parameters
    ----------
    mod : module
        An instance of a module.
    pattern : str, optional
        Only get functions whose name match a regular expression.

    Returns
    -------
    dict
        Dictionary of functions contained in the module. Keys are the
        function names, values are the functions themselves.
    """
    import inspect
    import re

    funcs = {}
    for name, func in inspect.getmembers(mod, inspect.isroutine):
        if pattern is None or re.match(pattern, name):
            funcs[name] = func
    return funcs


def add_functions_to_class(cls, funcs):
    """Add functions as methods of a class.

    Parameters
    ----------
    cls : class
        A class.
    funcs : dict
        Dictionary of function names and instances.
    """
    for name, func in funcs.items():
        setattr(cls, name, func)


def add_module_functions_to_class(cls, module, pattern=None):
    """Add functions from a module to a class as methods.

    Parameters
    ----------
    cls : class
        A class.
    module : module
        An instance of a module.
    pattern : str, optional
        Only get functions whose name match a regular expression.
    """
    import inspect
    import imp
    import os
    
    caller = inspect.stack()[1]
    path = os.path.join(os.path.dirname(caller[1]), os.path.dirname(module))

    (module, _) = os.path.splitext(os.path.basename(module))

    mod = imp.load_module(module, *imp.find_module(module, [path]))

    funcs = get_functions_from_module(mod, pattern=pattern)
    strip_grid_from_method_docstring(funcs)
    add_functions_to_class(cls, funcs)


def strip_grid_from_method_docstring(funcs):
    """Remove 'grid' from the parameters of a dict of functions' docstrings.

    Note that the docstring must be close to numpydoc standards for this to
    work.

    Parameters
    ----------
    funcs : dict
        Dictionary of functions to modify. Keys are the function names,
        values are the functions themselves.

    Examples
    --------
    >>> from landlab.grid.mappers import dummy_func_to_demonstrate_docstring_modification as dummy_func
    >>> funcs = {'dummy_func_to_demonstrate_docstring_modification':
    ...          dummy_func}
    >>> help(dummy_func)
    Help on function dummy_func_to_demonstrate_docstring_modification in module landlab.grid.mappers:
    <BLANKLINE>
    dummy_func_to_demonstrate_docstring_modification(grid, some_arg)
        A dummy function to demonstrate automated docstring changes.
    <BLANKLINE>
        Construction::
    <BLANKLINE>
            dummy_func_to_demonstrate_docstring_modification(grid, some_arg)
    <BLANKLINE>
        Parameters
        ----------
        grid : ModelGrid
            A Landlab modelgrid.
        some_arg : whatever
            A dummy argument.
    <BLANKLINE>
        Examples
        --------
        ...
    <BLANKLINE>
    >>> strip_grid_from_method_docstring(funcs)
    >>> help(dummy_func)
    Help on function dummy_func_to_demonstrate_docstring_modification in module landlab.grid.mappers:
    <BLANKLINE>
    dummy_func_to_demonstrate_docstring_modification(grid, some_arg)
        A dummy function to demonstrate automated docstring changes.
    <BLANKLINE>
        Construction::
    <BLANKLINE>
            grid.dummy_func_to_demonstrate_docstring_modification(some_arg)
    <BLANKLINE>
        Parameters
        ----------
        some_arg : whatever
            A dummy argument.
    <BLANKLINE>
        Examples
        --------
        ...
    <BLANKLINE>
    """
    import re
    for name, func in funcs.items():
        # strip the entry under "Parameters":
        func.__doc__ = re.sub('grid *:.*?\n.*?\n *', '', func.__doc__)
        # # cosmetic magic to get a two-line signature to line up right:
        match_2_lines = re.search(func.__name__+'\(grid,[^\)]*?\n.*?\)',
                                  func.__doc__)
        try:
            lines_were = match_2_lines.group()
        except AttributeError:  # no successful match
            pass
        else:
            end_chars = re.search('    .*?\)', lines_were).group()[4:]
            lines_are_now = re.sub('    .*?\)', '         '+end_chars,
                                   lines_were)
            func.__doc__ = (func.__doc__[:match_2_lines.start()] +
                            lines_are_now +
                            func.__doc__[match_2_lines.end():])
        # Move "grid" in signature from an arg to the class position
        func.__doc__ = re.sub(func.__name__+'\(grid, ',
                              'grid.'+func.__name__+'(',
                              func.__doc__)


def argsort_points_by_x_then_y(pts):
    """Sort points by coordinates, first x then y, returning sorted indices.

    Parameters
    ----------
    pts : Nx2 NumPy array of float
        (x,y) points to be sorted

    Returns
    -------
    indices : Nx1 NumPy array of int
        indices of sorted (x,y) points
    
    Examples
    --------
    >>> import numpy as np
    >>> from landlab.core.utils import argsort_points_by_x_then_y
    >>> pts = np.zeros((10, 2))
    >>> pts[:,0] = np.array([0., 0., 0., 1., 1., 1., 1., 2., 2., 2.])
    >>> pts[:,1] = np.array([0., 1., 2., -0.5, 0.5, 1.5, 2.5, 0., 1., 2.])
    >>> idx = argsort_points_by_x_then_y(pts)
    >>> idx
    array([3, 0, 7, 4, 1, 8, 5, 2, 9, 6])
    """
    a = pts[:,0].argsort(kind='mergesort')
    b = pts[a,1].argsort(kind='mergesort')
    return as_id_array(a[b])


def sort_points_by_x_then_y(pts):
    """Sort points by coordinates, first x then y.

    Parameters
    ----------
    pts : Nx2 NumPy array of float
        (x,y) points to be sorted

    Returns
    -------
    pts : Nx2 NumPy array of float
        sorted (x,y) points

    Examples
    --------
    >>> import numpy as np
    >>> from landlab.core.utils import sort_points_by_x_then_y
    >>> pts = np.zeros((10, 2))
    >>> pts[:,0] = np.array([0., 0., 0., 1., 1., 1., 1., 2., 2., 2.])
    >>> pts[:,1] = np.array([0., 1., 2., -0.5, 0.5, 1.5, 2.5, 0., 1., 2.])
    >>> pts = sort_points_by_x_then_y(pts)
    >>> pts
    array([[ 1. , -0.5],
           [ 0. ,  0. ],
           [ 2. ,  0. ],
           [ 1. ,  0.5],
           [ 0. ,  1. ],
           [ 2. ,  1. ],
           [ 1. ,  1.5],
           [ 0. ,  2. ],
           [ 2. ,  2. ],
           [ 1. ,  2.5]])
    """
    indices = argsort_points_by_x_then_y(pts)
    pts[:, 0] = pts[indices, 0]
    pts[:, 1] = pts[indices, 1]
    return pts


def anticlockwise_argsort_points(pts, midpt=None):
    """Argort points into anticlockwise order around a supplied center.

    Sorts CCW from east. Assumes a convex hull.

    Parameters
    ----------
    pts : Nx2 NumPy array of float
        (x,y) points to be sorted
    midpt : len-2 NumPy array of float (optional)
        (x, y) of point about which to sort. If not provided, mean of pts is
        used.

    Returns
    -------
    pts : Nx2 NumPy array of float
        sorted (x,y) points

    Examples
    --------
    >>> import numpy as np
    >>> from landlab.core.utils import anticlockwise_argsort_points
    >>> pts = np.zeros((4, 2))
    >>> pts[:,0] = np.array([-3., -1., -1., -3.])
    >>> pts[:,1] = np.array([-1., -3., -1., -3.])
    >>> sortorder = anticlockwise_argsort_points(pts)
    >>> np.all(sortorder == np.array([2, 0, 3, 1]))
    True
    """
    if midpt is None:
        midpt = pts.mean(axis=0)
    assert len(midpt) == 2
    theta = np.arctan2(pts[:, 1] - midpt[1], pts[:, 0] - midpt[0])
    theta = theta % (2.*np.pi)
    sortorder = np.argsort(theta)
    return sortorder


if __name__ == '__main__':
    import doctest
    doctest.testmod()

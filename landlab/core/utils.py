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
    ~landlab.core.utils.get_categories_from_grid_methods
"""
from __future__ import print_function

import numpy as np

SIZEOF_INT = np.dtype(np.int).itemsize


def degrees_to_radians(degrees):
    """Convert compass-style degrees to radians.

    Convert angles in degrees measured clockwise starting from north to
    angles measured counter-clockwise from the positive x-axis in radians

    Parameters
    ----------
    degrees : float or ndarray
        Converted angles in degrees.

    Returns
    -------
    rads : float or ndarray
        Angles in radians.

    Examples
    --------
    >>> import numpy as np
    >>> from landlab.core.utils import degrees_to_radians

    >>> degrees_to_radians(90.0)
    0.0
    >>> degrees_to_radians(0.0) == np.pi / 2.
    True
    >>> degrees_to_radians(-180.0) == 3. * np.pi / 2.
    True
    >>> np.testing.assert_array_almost_equal([ np.pi, np.pi],
    ...                                       degrees_to_radians([ -90.,  270.]))
    """
    rads = np.pi * np.array(degrees) / 180.0

    return (5.0 * np.pi / 2.0 - rads) % (2.0 * np.pi)


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
    degrees = (5.0 * np.pi / 2.0 - rads) % (2.0 * np.pi)
    return 180.0 / np.pi * degrees


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
    if hasattr(x, "ext"):
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

    >>> x = np.arange(5, dtype=np.intp)
    >>> y = np.where(x < 3)[0]
    >>> y.dtype == np.intp
    True
    >>> as_id_array(y).dtype == np.int
    True
    """
    try:
        if array.dtype == np.int:
            return array.view(np.int)
        else:
            return array.astype(np.int)
    except AttributeError:
        return np.asarray(array, dtype=np.int)


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
        ids = as_id_array(np.asarray(args[0])).reshape((-1,))
    else:
        raise ValueError("Number of arguments must be 0 or 1.")

    return ids


def get_functions_from_module(mod, pattern=None, exclude=None):
    """Get all the function in a module.

    Parameters
    ----------
    mod : module
        An instance of a module.
    pattern : str, optional
        Only get functions whose name match a regular expression.
    exclude : str, optional
        Only get functions whose name exclude the regular expression.

    *Note* if both pattern and exclude are provided both conditions must be met.

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
        if pattern is None or re.search(pattern, name):
            if exclude is None or (re.search(exclude, name) is None):
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


def add_module_functions_to_class(cls, module, pattern=None, exclude=None):
    """Add functions from a module to a class as methods.

    Parameters
    ----------
    cls : class
        A class.
    module : module
        An instance of a module.
    pattern : str, optional
        Only get functions whose name match a regular expression.
    exclude : str, optional
        Only get functions whose name exclude the regular expression.

    *Note* if both pattern and exclude are provided both conditions must be met.
    """
    import inspect
    import imp
    import os

    caller = inspect.stack()[1]
    path = os.path.join(os.path.dirname(caller[1]), os.path.dirname(module))

    (module, _) = os.path.splitext(os.path.basename(module))

    mod = imp.load_module(module, *imp.find_module(module, [path]))

    funcs = get_functions_from_module(mod, pattern=pattern, exclude=exclude)
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
        func.__doc__ = re.sub("grid *:.*?\n.*?\n *", "", func.__doc__)
        # # cosmetic magic to get a two-line signature to line up right:
        match_2_lines = re.search(
            func.__name__ + r"\(grid,[^\)]*?\n.*?\)", func.__doc__
        )
        try:
            lines_were = match_2_lines.group()
        except AttributeError:  # no successful match
            pass
        else:
            end_chars = re.search(r"    .*?\)", lines_were).group()[4:]
            lines_are_now = re.sub(r"    .*?\)", "         " + end_chars, lines_were)
            func.__doc__ = (
                func.__doc__[: match_2_lines.start()]
                + lines_are_now
                + func.__doc__[match_2_lines.end() :]
            )
        # Move "grid" in signature from an arg to the class position
        func.__doc__ = re.sub(
            func.__name__ + r"\(grid, ", "grid." + func.__name__ + "(", func.__doc__
        )


def argsort_points_by_x_then_y(points):
    """Sort points by coordinates, first x then y, returning sorted indices.

    Parameters
    ----------
    points : tuple of ndarray or ndarray of float, shape `(*, 2)`
        Coordinates of points to be sorted. Sort by first coordinate, then
        second.

    Returns
    -------
    ndarray of int, shape `(n_points, )`
        Indices of sorted points.

    Examples
    --------
    >>> import numpy as np
    >>> from landlab.core.utils import argsort_points_by_x_then_y

    >>> points = np.zeros((10, 2))
    >>> points[:, 0] = np.array([0., 0., 0.,  1.,  1.,  1.,  1.,  2., 2., 2.])
    >>> points[:, 1] = np.array([0., 1., 2., -0.5, 0.5, 1.5, 2.5, 0., 1., 2.])
    >>> argsort_points_by_x_then_y(points)
    array([3, 0, 7, 4, 1, 8, 5, 2, 9, 6])

    >>> x = [0., 0., 0.,
    ...      1., 1., 1., 1.,
    ...      2., 2., 2.]
    >>> y = [ 0. , 1. , 2. ,
    ...      -0.5, 0.5, 1.5, 2.5,
    ...       0. , 1. , 2.]
    >>> indices = argsort_points_by_x_then_y((x, y))
    >>> indices
    array([3, 0, 7, 4, 1, 8, 5, 2, 9, 6])

    >>> argsort_points_by_x_then_y(np.array((x, y)))
    array([3, 0, 7, 4, 1, 8, 5, 2, 9, 6])
    """
    if isinstance(points, np.ndarray):
        if points.shape[0] > points.shape[1]:
            points = points.T
        try:
            return argsort_points_by_x_then_y((points[0, :], points[1, :]))
        except IndexError:
            return as_id_array([0])
    else:
        points = [np.asarray(coord) for coord in points]
        a = points[0].argsort(kind="mergesort")
        b = points[1][a].argsort(kind="mergesort")
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
        pts : N NumPy array of int
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
    theta = theta % (2.0 * np.pi)
    sortorder = np.argsort(theta)
    return sortorder


def anticlockwise_argsort_points_multiline(pts_x, pts_y, out=None):
    """Argort multi lines of points into CCW order around the geometric center.

    This version sorts columns of data in a 2d array. Sorts CCW from east
    around the geometric center of the points in the row.
    Assumes a convex hull.

    Parameters
    ----------
    pts_x : rows x n_elements array of float
        rows x points_to_sort x x_coord of points
    pts_y : rows x n_elements array of float
        rows x points_to_sort x y_coord of points
    out : rows x n_elements (optional)
        If provided, the ID array to be sorted

    Returns
    -------
    sortorder : rows x n_elements NumPy array of int
        sorted (x,y) points

    Examples
    --------
    >>> import numpy as np
    >>> from landlab.core.utils import anticlockwise_argsort_points_multiline
    >>> pts = np.array([[1, 3, 0, 2], [2, 0, 3, 1]])
    >>> pts_x = np.array([[-3., -1., -1., -3.], [-3., -1., -1., -3.]])
    >>> pts_y = np.array([[-1., -3., -1., -3.], [-3., -1., -3., -1.]])
    >>> sortorder = anticlockwise_argsort_points_multiline(
    ...     pts_x, pts_y, out=pts)
    >>> np.all(sortorder == np.array([[2, 0, 3, 1], [1, 3, 0, 2]]))
    True
    >>> np.all(pts == np.array([[0, 1, 2, 3], [0, 1, 2, 3]]))
    True
    """
    nrows = pts_x.shape[0]
    midpt = np.empty((nrows, 2), dtype=float)
    midpt[:, 0] = pts_x.mean(axis=1)
    midpt[:, 1] = pts_y.mean(axis=1)
    theta = np.arctan2(
        pts_y - midpt[:, 1].reshape((nrows, 1)), pts_x - midpt[:, 0].reshape((nrows, 1))
    )
    theta = theta % (2.0 * np.pi)
    sortorder = np.argsort(theta)
    if out is not None:
        out[:] = out[np.ogrid[:nrows].reshape((nrows, 1)), sortorder]
    return sortorder


def get_categories_from_grid_methods(grid_type):
    """
    Create a dict of category:[method_names] for a LL grid type.

    Looks in the final line of the docstrings
    of class methods and properties for a catgory declaration, "LLCATS: ".
    It then creates and returns a dict with keys found as categories and
    values that are lists of the names of methods that have that category.

    Currently defined LLCATS are:
        DEPR : deprecated
        GINF : information about the grid as a whole
        NINF : information about nodes
        LINF : information about links
        PINF : information about patches
        CINF : information about cells
        FINF : information about faces
        CNINF : information about corners
        FIELDIO : methods to access and change fields
        FIELDADD : methods to create new fields/delete old fields
        FIELDINF : information about fields (keys, names, etc)
        GRAD : methods for gradients, fluxes, divergences and slopes
        MAP : methods to map from one element type to another
        BC : methods to interact with BCs
        SURF : methods for surface analysis (slope, aspect, hillshade)
        SUBSET : methods to indentify part of the grid based on conditions
        CONN : method describing the connectivity of one element to another
               (i.e., 'links_at_node')
        MEAS : method describing a quantity defined on an element (i.e.,
               'length_of_link')
        OTHER : anything else

    Parameters
    ----------
    grid_type : {'ModelGrid', 'RasterModelGrid', 'HexModelGrid',
                 'RadialModelGrid', 'VoronoiDelaunayGrid', 'NetworkModelGrid'}
        String of raster to inspect.

    Returns
    -------
    cat_dict : dict
        Dictionary with cats as keys and lists of method name strings as
        values.
    grid_dict : dict
        Dictionary with method name strings as keys and lists of cats as
        values.
    FAILS : dict of dicts
        contains any problematic LLCAT entries. Keys: 'MISSING' - list of names
        of any public method or property without an LLCAT declared.
    """
    import inspect
    import re
    from landlab import (
        ModelGrid,
        RasterModelGrid,
        HexModelGrid,
        RadialModelGrid,
        VoronoiDelaunayGrid,
        NetworkModelGrid,
    )
    from copy import copy

    grid_str_to_grid = {
        "ModelGrid": ModelGrid,
        "RasterModelGrid": RasterModelGrid,
        "HexModelGrid": HexModelGrid,
        "RadialModelGrid": RadialModelGrid,
        "VoronoiDelaunayGrid": VoronoiDelaunayGrid,
        "NetworkModelGrid": NetworkModelGrid,
    }
    grid_dict = {}
    cat_dict = {}
    FAILS = {"MISSING": []}
    grid = grid_str_to_grid[grid_type]
    funcs = {}
    for name, func in inspect.getmembers(grid):
        funcs[name] = func
    for method_name in funcs.keys():
        if method_name[0] == "_":
            continue
        else:
            method_doc = funcs[method_name].__doc__
            try:
                cat_str = re.search("LLCATS:.+", method_doc)
            except TypeError:
                pass
            else:
                if cat_str is None:
                    FAILS["MISSING"].append(method_name)
                    continue
                cats = cat_str.group().split()[1:]
                grid_dict[method_name] = copy(cats)
                for cat in cats:
                    try:
                        cat_dict[cat].append(method_name)
                    except KeyError:
                        cat_dict[cat] = [method_name]

    return cat_dict, grid_dict, FAILS


if __name__ == "__main__":
    import doctest

    doctest.testmod()

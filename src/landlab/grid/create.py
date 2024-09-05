#! /usr/bin/env python
"""Create landlab model grids."""

from ..core import load_params
from ..io import read_esri_ascii
from ..io.netcdf import read_netcdf
from ..values import constant
from ..values import plane
from ..values import random
from ..values import sine
from ..values import units
from .hex import HexModelGrid
from .network import NetworkModelGrid
from .radial import RadialModelGrid
from .raster import RasterModelGrid
from .voronoi import VoronoiDelaunayGrid

_MODEL_GRIDS = {
    "RasterModelGrid": RasterModelGrid,
    "HexModelGrid": HexModelGrid,
    "VoronoiDelaunayGrid": VoronoiDelaunayGrid,
    "NetworkModelGrid": NetworkModelGrid,
    "RadialModelGrid": RadialModelGrid,
}

_SYNTHETIC_FIELD_CONSTRUCTORS = {
    "plane": plane,
    "random": random,
    "sine": sine,
    "constant": constant,
    "units": units,
}


class Error(Exception):
    """Base class for exceptions from this module."""

    pass


class BadGridTypeError(Error):
    """Raise this error for a bad grid type."""

    def __init__(self, grid_type):
        self._type = str(grid_type)  # TODO: not tested.

    def __str__(self):
        return self._type  # TODO: not tested.


def grid_from_dict(grid_type, params):
    """Create a grid from a dictionary of parameters."""
    try:
        cls = _MODEL_GRIDS[grid_type]
    except KeyError as exc:
        raise ValueError(f"unknown grid type ({grid_type})") from exc
    args, kwargs = _parse_args_kwargs(params)
    return cls(*args, **kwargs)


def grids_from_file(file_like, section=None):
    """Create grids from a file."""
    params = load_params(file_like)

    if section:
        try:
            grids = params[section]
        except KeyError as exc:  # TODO: not tested.
            raise ValueError(
                f"missing required section ({section})"
            ) from exc  # TODO: not tested.
    else:  # TODO: not tested.
        grids = params  # TODO: not tested.

    new_grids = []
    for grid_type, grid_desc in as_list_of_tuples(grids):
        new_grids.append(grid_from_dict(grid_type, grid_desc))

    return new_grids


def add_fields_from_dict(grid, fields):
    """Add fields to a grid from a dictionary."""
    fields = dict(fields)

    unknown_locations = set(fields) - set(grid.VALID_LOCATIONS)
    if unknown_locations:
        raise ValueError(
            "unknown field locations ({})".format(", ".join(unknown_locations))
        )

    for location, fields_at_location in fields.items():
        for name, function in fields_at_location.items():
            add_field_from_function(grid, name, function, at=location)

    return grid


def add_field_from_function(grid, name, functions, at="node"):
    """Add a field to a grid as functions.

    Parameters
    ----------
    grid : ModelGrid
        A landlab grid to add fields to.
    name : str
        Name of the new field.
    functions : *(func_name, func_args)* or iterable of *(func_name, func_args)*
        The functions to apply to the field. Functions are applied in the order
        the appear in the list.
    at : str
        The grid element to which the field will be added.

    Returns
    -------
    ModelGrid
        The grid with the new field.
    """
    valid_functions = set(_SYNTHETIC_FIELD_CONSTRUCTORS) | {
        "read_esri_ascii",
        "read_netcdf",
    }

    for func_name, func_args in as_list_of_tuples(functions):
        if func_name not in valid_functions:
            raise ValueError(f"function not understood ({func_name})")

        args, kwargs = _parse_args_kwargs(func_args)

        if func_name in _SYNTHETIC_FIELD_CONSTRUCTORS:
            # if any args, raise an error, there shouldn't be any.
            synth_function = _SYNTHETIC_FIELD_CONSTRUCTORS[func_name]
            synth_function(grid, name, at=at, **kwargs)
        elif func_name == "read_esri_ascii":
            read_esri_ascii(*args, grid=grid, name=name, **kwargs)
        elif func_name == "read_netcdf":
            read_netcdf(*args, grid=grid, name=name, **kwargs)

    return grid


def add_boundary_conditions(grid, boundary_conditions=()):
    for bc_name, bc_args in as_list_of_tuples(boundary_conditions):
        args, kwargs = _parse_args_kwargs(bc_args)
        try:
            func = getattr(grid, bc_name)
        except AttributeError as exc:
            raise ValueError(
                f"create_grid: No function {bc_name} exists for grid types "
                f"{grid.__class__.__name__}. If you think this type of grid "
                "should have such a function. Please create a GitHub Issue to "
                "discuss contributing it to the Landlab codebase."
            ) from exc
        else:
            func(*args, **kwargs)


def as_list_of_tuples(items):
    """Convert a collection of key/values to a list of tuples.

    Examples
    --------
    >>> from collections import OrderedDict
    >>> from landlab.grid.create import as_list_of_tuples
    >>> as_list_of_tuples({"eric": "idle"})
    [('eric', 'idle')]
    >>> as_list_of_tuples([("john", "cleese"), {"eric": "idle"}])
    [('john', 'cleese'), ('eric', 'idle')]
    >>> as_list_of_tuples(
    ...     [("john", "cleese"), OrderedDict([("eric", "idle"), ("terry", "gilliam")])]
    ... )
    [('john', 'cleese'), ('eric', 'idle'), ('terry', 'gilliam')]
    """
    try:
        items = list(items.items())
    except AttributeError:
        items = list(items)

    if len(items) == 2 and isinstance(items[0], str):
        items = [items]

    tuples = []
    for item in items:
        try:
            tuples.extend(list(item.items()))
        except AttributeError:
            tuples.append(tuple(item))

    return tuples


def create_grid(file_like, section=None):
    """Create grid, initialize fields, and set boundary conditions.

    **create_grid** expects a dictionary with three keys: "grid", "fields", and
    "boundary_conditions".

    **Dictionary Section "grid"**

    The value associated with the "grid" key should itself be a dictionary
    containing the name of a Landlab model grid type as its only key. The
    following grid types are valid:

        -  :py:class:`~landlab.grid.raster.RasterModelGrid`
        -  :py:class:`~landlab.grid.voronoi.VoronoiDelaunayGrid`
        -  :py:class:`~landlab.grid.hex.HexModelGrid`
        -  :py:class:`~landlab.grid.radial.RadialModelGrid`
        -  :py:class:`~landlab.grid.network.NetworkModelGrid`

    The value associated with the grid name key is a list containing the
    arguments. If any keyword arguments are passed, they should be passed as
    the last element of the list. For example the following code block is a
    yaml file indicating a RasterModelGrid with shape (4, 5) and xy-spacing of
    (3, 4).

    .. code-block:: yaml

        grid:
          RasterModelGrid:
            - [4, 5]
            - xy_spacing: [3, 4]

    These arguments and keyword arguments will be passed to the ``__init__``
    constructor of the specified model grid. Refer to the documentation for
    each grid to determine its requirements.

    **Dictionary Section "fields"**

    Fields can be created by reading from files or by creating synthetic
    values.

    The value associated with the "fields" key is a nested set of dictionaries
    indicating where the fields are created, what the field names are, and how
    to create the fields. As part of a grid's description, the value
    associated with the "fields" key must be a dictionary with keys indicating
    at which grid elements fields should be created (e.g. to create fields at
    node, use "node").

    The value associated with each "xxx" (i.e. "node", "link", "patch", etc.)
    value is itself a dictionary
    indicating the name of the field and how it should be created. A field can
    either be created by reading from a file or creating synthetic values. The
    :py:func:`~landlab.io.netcdf.read.read_netcdf` and
    :py:func:`~landlab.io.esri_ascii.read_esri_ascii` functions, and the
    :py:mod:`synthetic fields <landlab.values.synthetic>`
    package are currently supported methods to create fields. These may be
    chained together (as is shown in the Example section below). If these
    functions do not meet your needs, we welcome contributions that extend the
    capabilities of this function.

    An additional supported method, which can be chained together with
    either synthetic fields or fields read from a file, is *units*. The
    *units* function will set the units attribute of its corresponding field.
    If this optional function is not used, the resulting field will not be
    given any units.

    The following example would use the
    :py:func:`~landlab.values.synthetic.plane` function from the synthetic
    values package to create a *node* value for the field
    *topographic__elevation*. The plane function adds values to a Landlab model
    grid field that lie on a plane specified by a point and a normal vector. In
    the below example the plane goes through the point (1.0, 1.0, 1.0) and has
    a normal of (-2.0, -1.0, 1.0). The *units* function sets the units of
    elevation to *meters*.

    .. code-block:: yaml

        grid:
          RasterModelGrid:
            - [4, 5]
            - xy_spacing: [3, 4]
            - fields:
                node:
                  topographic__elevation:
                    plane:
                      - point: [1, 1, 1]
                        normal: [-2, -1, 1]
                    units:
                        - units: "meters"

    **Dictionary Section "boundary_conditions"**

    The final portion of the input dictionary calls bound functions of the
    model grid to set boundary conditions. Any valid bound function can be
    called. The specified functions are provided in a list, and called in
    order. If required, multiple functions may be called.

    Each entry to the list is a dictionary with a single key, the name of the
    bound function. The value associated with that key is a list of arguments
    and keyword arguments, similar in structure to those described above.
    As with the "fields" section, the "boundary_conditions" section must be
    described under its associated grid description.

    For example, the following sets closed boundaries at all sides of the grid.

    .. code-block:: yaml

        grid:
          RasterModelGrid:
            - [4, 5]
            - xy_spacing: [3, 4]
            - boundary_conditions:
              - set_closed_boundaries_at_grid_edges:
                - True
                - True
                - True
                - True

    Parameters
    ----------
    file_like : file_like or str
        Dictionary, contents of a dictionary as a string, a file-like object,
        or the path to a file containing a YAML dictionary.

    Examples
    --------
    >>> import numpy as np
    >>> from landlab import create_grid
    >>> np.random.seed(42)
    >>> p = {
    ...     "grid": {
    ...         "RasterModelGrid": [
    ...             (4, 5),
    ...             {"xy_spacing": (3, 4)},
    ...             {
    ...                 "fields": {
    ...                     "node": {
    ...                         "spam": {
    ...                             "plane": [
    ...                                 {"point": (1, 1, 1), "normal": (-2, -1, 1)}
    ...                             ],
    ...                             "random": [
    ...                                 {"distribution": "uniform", "low": 1, "high": 4}
    ...                             ],
    ...                         }
    ...                     },
    ...                     "link": {
    ...                         "eggs": {
    ...                             "constant": [{"where": "ACTIVE_LINK", "value": 12}]
    ...                         }
    ...                     },
    ...                 }
    ...             },
    ...             {
    ...                 "boundary_conditions": [
    ...                     {
    ...                         "set_closed_boundaries_at_grid_edges": [
    ...                             True,
    ...                             True,
    ...                             True,
    ...                             True,
    ...                         ]
    ...                     }
    ...                 ]
    ...             },
    ...         ]
    ...     }
    ... }
    >>> mg = create_grid(p, section="grid")
    >>> mg.number_of_nodes
    20
    >>> "spam" in mg.at_node
    True
    >>> "eggs" in mg.at_link
    True
    >>> mg.x_of_node
    array([  0.,   3.,   6.,   9.,  12.,
             0.,   3.,   6.,   9.,  12.,
             0.,   3.,   6.,   9.,  12.,
             0.,   3.,   6.,   9.,  12.])
    >>> mg.status_at_node
    array([4, 4, 4, 4, 4,
           4, 0, 0, 0, 4,
           4, 0, 0, 0, 4,
           4, 4, 4, 4, 4], dtype=uint8)
    >>> np.round(mg.at_node["spam"].reshape(mg.shape), decimals=2)
    array([[ 0.12,   7.85,  13.2 ,  18.8 ,  23.47],
           [ 3.47,   9.17,  17.6 ,  22.8 ,  29.12],
           [ 7.06,  15.91,  21.5 ,  25.64,  31.55],
           [11.55,  17.91,  24.57,  30.3 ,  35.87]])
    """
    if isinstance(file_like, dict):
        params = file_like
    else:
        params = load_params(file_like)

    if section:
        grids = params[section]
    else:
        grids = params

    new_grids = []
    for grid_type, grid_desc in as_list_of_tuples(grids):
        grid_desc = norm_grid_description(grid_desc)

        fields = grid_desc.pop("fields", {})
        boundary_conditions = grid_desc.pop("boundary_conditions", {})

        grid = grid_from_dict(grid_type, grid_desc)
        add_fields_from_dict(grid, fields)
        add_boundary_conditions(grid, boundary_conditions)

        new_grids.append(grid)

    if len(new_grids) == 1:
        return new_grids[0]
    else:
        return new_grids


def norm_grid_description(grid_desc):
    """Normalize a grid description into a canonical form.

    Examples
    --------
    >>> from landlab.grid.create import norm_grid_description

    >>> grid_desc = [(3, 4), {"xy_spacing": 4.0, "xy_of_lower_left": (1.0, 2.0)}]
    >>> normed_items = list(norm_grid_description(grid_desc).items())
    >>> normed_items.sort()
    >>> normed_items
    [('args', [(3, 4)]), ('xy_of_lower_left', (1.0, 2.0)), ('xy_spacing', 4.0)]
    """
    if not isinstance(grid_desc, dict):
        args, kwds = [], {}
        for arg in grid_desc:
            if isinstance(arg, dict) and {"fields", "boundary_conditions"} & set(
                arg.keys()
            ):
                kwds.update(arg)
            else:
                args.append(arg)
        if isinstance(args[-1], dict):
            kwds.update(args.pop())
        kwds.update({"args": args})
        return kwds
    return grid_desc


def _parse_args_kwargs(list_of_args_kwargs):
    if isinstance(list_of_args_kwargs, dict):
        args, kwargs = list_of_args_kwargs.pop("args", ()), list_of_args_kwargs
        if not isinstance(args, (tuple, list)):
            args = (args,)
    else:
        args, kwargs = [], {}
        for arg in list(list_of_args_kwargs):
            if isinstance(arg, dict) and {"fields", "boundary_conditions"} & set(
                arg.keys()
            ):
                kwargs.update(arg)  # TODO: not tested.
            else:
                args.append(arg)
        if isinstance(args[-1], dict):
            kwargs.update(args.pop())

    return tuple(args), kwargs

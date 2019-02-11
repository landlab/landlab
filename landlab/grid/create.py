#! /usr/bin/env python
"""Create landlab model grids."""
import inspect
from warnings import warn

from landlab.core import load_params, model_parameter_dictionary as mpd
from landlab.io import read_esri_ascii
from landlab.io.netcdf import read_netcdf
from landlab.values import constant, plane, random, sine

from .hex import HexModelGrid, from_dict as hex_from_dict
from .network import NetworkModelGrid
from .radial import RadialModelGrid
from .raster import RasterModelGrid, from_dict as raster_from_dict
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
}


class Error(Exception):

    """Base class for exceptions from this module."""

    pass


class BadGridTypeError(Error):

    """Raise this error for a bad grid type."""

    def __init__(self, grid_type):
        self._type = str(grid_type)

    def __str__(self):
        return self._type


_GRID_READERS = {"raster": raster_from_dict, "hex": hex_from_dict}


def create_and_initialize_grid(input_source):
    """Create and initialize a grid from a file.

    Creates, initializes, and returns a new grid object using parameters
    specified in *input_source*. *input_source* is either a
    ModelParameterDictionary instance (or, really, just dict-like) or a
    named input file.

    Parameters
    ----------
    input_source : str or dict
        Input file or ``dict`` of parameter values.

    Raises
    ------
    KeyError
        If missing a key from the input file.

    Examples
    --------
    >>> from six import StringIO
    >>> import pytest
    >>> test_file = StringIO('''
    ... GRID_TYPE:
    ... raster
    ... NUM_ROWS:
    ... 4
    ... NUM_COLS:
    ... 5
    ... GRID_SPACING:
    ... 2.5
    ... ''')
    >>> from landlab import create_and_initialize_grid
    >>> with pytest.deprecated_call():
    ...    grid = create_and_initialize_grid(test_file)
    >>> grid.number_of_nodes
    20
    """
    msg = (
        "create_and_initialize_grid is deprecated and will be removed "
        "in landlab 2.0. Use create_grid instead."
    )
    warn(msg, DeprecationWarning)
    if isinstance(input_source, dict):
        param_dict = input_source
    else:
        param_dict = mpd.ModelParameterDictionary(from_file=input_source)

    grid_type = param_dict["GRID_TYPE"]

    grid_type.strip().lower()

    # Read parameters appropriate to that type, create it, and initialize it
    try:
        grid_reader = _GRID_READERS[grid_type]
    except KeyError:
        raise BadGridTypeError(grid_type)

    # Return the created and initialized grid
    return grid_reader(param_dict)


def create_grid(file_like):
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
    to create the fields. At the highest hierachical level, the value
    associated with the "fields" key must be a dictionary with keys indicating
    at which grid elements fields should be created (e.g. to create fields at
    node, use "at_node").

    The value associated with each "at_xxx" value is itself a dictionary
    indicating the name of the field an how it should be created. A field can
    either be created by reading from a file or creating synthetic values. The
    :py:func:`~landlab.io.netcdf.read.read_netcdf` and
    :py:func:`~landlab.io.esri_ascii.read_esri_ascii` functions, and the
    :py:mod:`synthetic fields <landlab.values.synthetic>`
    package are currently supported methods to create fields. These may be
    chained together (as is shown in the Example section below). If these
    functions do not meet your needs, we welcome contributions that extend the
    capabilities of this function.

    The following example would uses the
    :py:func:`~landlab.values.synthetic.plane` function from the synthetic
    values package to create an at_node value for the field
    topographic__elevation. The plane function adds values to a Landlab model
    grid field that lie on a plane specified by a point and a normal vector. In
    the below example the plane goes through the point (1.0, 1.0, 1.0) and has
    a normal of (-2.0, -1.0, 1.0).

    .. code-block:: yaml

        fields:
          at_node:
            topographic__elevation:
              plane:
                - point: [1, 1, 1]
                  normal: [-2, -1, 1]

    **Dictionary Section "boundary_conditions"**

    The final portion of the input dictionary calls bound functions of the
    model grid to set boundary conditions. Any valid bound function can be
    called. The specified functions are provided in a list, and called in
    order. If required, multiple functions may be called.

    Each entry to the list is a dictionary with a single key, the name of the
    bound function. The value associated with that key is a list of arguments
    and keyword arguments, similar in structure to those described above.

    For example, the following sets closed boundaries at all sides of the grid.

    .. code-block:: yaml

        boundary_conditions:
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
    >>> np.random.seed(42)
    >>> from landlab import create_grid
    >>> p = {'grid': {'RasterModelGrid': [(4,5),
    ...                                   {'xy_spacing': (3, 4)}]
    ...               },
    ...      'fields': {'at_node': {'spam': {'plane': [{'point': (1, 1, 1),
    ...                                                'normal': (-2, -1, 1)}],
    ...                                      'random': [{'distribution': 'uniform',
    ...                                                 'low': 1,
    ...                                                 'high': 4}]
    ...                                      },
    ...                             },
    ...                 'at_link': {'eggs': {'constant': [{'where': 'ACTIVE_LINK',
    ...                                                    'constant': 12}],
    ...                                       },
    ...                             },
    ...                 },
    ...      'boundary_conditions': [
    ...                     {'set_closed_boundaries_at_grid_edges':
    ...                                         [True, True, True, True]
    ...                      }]
    ...      }
    >>> mg = create_grid(p)
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
    >>> np.round(mg.at_node['spam'].reshape(mg.shape), decimals=2)
    array([[  0.12,   7.85,  13.2 ,  18.8 ,  23.47],
           [  3.47,   9.17,  17.6 ,  22.8 ,  29.12],
           [  7.06,  15.91,  21.5 ,  25.64,  31.55],
           [ 11.55,  17.91,  24.57,  30.3 ,  35.87]])
    """
    # part 0, parse input
    if isinstance(file_like, dict):
        dict_like = file_like
    else:
        dict_like = load_params(file_like)

    # part 1 create grid
    grid_dict = dict_like.pop("grid", None)
    if grid_dict is None:
        msg = "create_grid: no grid dictionary provided. This is required."
        raise ValueError(msg)

    for grid_type in grid_dict:
        if grid_type in _MODEL_GRIDS:
            grid_class = _MODEL_GRIDS[grid_type]
        else:
            msg = "create_grid: provided grid type not supported."
            raise ValueError

    if len(grid_dict) != 1:
        msg = (
            "create_grid: two entries to grid dictionary provided. "
            "This is not supported."
        )
        raise ValueError

    args, kwargs = _parse_args_kwargs(grid_dict.pop(grid_type))
    grid = grid_class(*args, **kwargs)

    # part two, create fields
    fields_dict = dict_like.pop("fields", {})

    # for each grid element:
    for at_group in fields_dict:
        at = at_group[3:]
        if at not in grid.groups:
            msg = (
                "create_grid: No field location ",
                "{at} ".format(at=at),
                "exists for grid types",
                "{grid}. ".format(grid=grid_type),
            )
            raise ValueError(msg)

        at_dict = fields_dict[at_group]
        # for field at grid element
        for name in at_dict:
            name_dict = at_dict[name]

            # for each function, add values.
            for func in name_dict:
                args, kwargs = _parse_args_kwargs(name_dict[func])
                if func in _SYNTHETIC_FIELD_CONSTRUCTORS:
                    # if any args, raise an error, there shouldn't be any.
                    synth_function = _SYNTHETIC_FIELD_CONSTRUCTORS[func]
                    synth_function(grid, name, at=at, **kwargs)
                elif func == "read_esri_ascii":
                    read_esri_ascii(*args, grid=grid, name=name, **kwargs)
                elif func == "read_netcdf":
                    read_netcdf(*args, grid=grid, name=name, **kwargs)
                else:
                    msg = (
                        "create_grid: Bad function ",
                        "{func} ".format(func=func),
                        "for creating a field.",
                    )
                    raise ValueError(msg)

    # part three, set boundary conditions
    bc_list = dict_like.pop("boundary_conditions", [])
    for bc_function_dict in bc_list:
        if len(bc_function_dict) != 1:
            msg = (
                "create_grid: two entries to a boundary condition function "
                "dictionary were provided. This is not supported."
            )
            raise ValueError(msg)
        for bc_function in bc_function_dict:
            args, kwargs = _parse_args_kwargs(bc_function_dict[bc_function])
            methods = dict(inspect.getmembers(grid, inspect.ismethod))
            if bc_function in methods:
                methods[bc_function](*args, **kwargs)
            else:
                msg = (
                    "create_grid: No function ",
                    "{func} ".format(func=bc_function),
                    "exists for grid types ",
                    "{grid}. ".format(grid=grid_type),
                    "If you think this type of grid should have such a ",
                    "function. Please create a GitHub Issue to discuss ",
                    "contributing it to the Landlab codebase.",
                )
                raise ValueError(msg)

    return grid


def _parse_args_kwargs(list_of_args_kwargs):
    if isinstance(list_of_args_kwargs[-1], dict):
        kwargs = list_of_args_kwargs.pop()
    else:
        kwargs = {}
    return list_of_args_kwargs, kwargs

#! /usr/bin/env python
"""Create landlab model grids."""
from warnings import warn

from landlab.core import model_parameter_dictionary as mpd
from landlab.io import read_esri_ascii, read_shapefile
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
    >>> with pytest.warns(DeprecationWarning):
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


def create_grid(dict_like):
    """
    """
    # part 1 create grid
    grid_dict = dict_like.pop("grid", None)
    if grid_dict is None:
        msg = "create_grid: no grid dictionary provided. This is required."
        raise ValueError(msg)

    grid_type = grid_dict.keys()[0]
    if grid_type in _MODEL_GRIDS:
        grid_class = _MODEL_GRIDS[grid_type]
    else:
        msg = "create_grid: provided grid type not supported."
        raise ValueError
    grid_params = grid_dict.pop(grid_type)
    if len(grid) is not 0:
        msg = (
            "create_grid: two entries to grid dictionary provided. "
            "This is not supported."
        )
        raise ValueError
    grid = grid_class.from_dict(grid_params)

    # part two, create fields
    fields_dict = dict_like.pop("fields", {})

    # for each grid element:
    for at in grid.groups:
        at_group = "at_" + grid
        at_dict = fields_dict.pop(at_group, None)

        # for field at grid element
        for name in at_dict:
            name_dict = at_dict.pop(name)

            # for each function, add values.
            for func in name_dict:
                if func in _SYNTHETIC_FIELD_CONSTRUCTORS:
                    synth_function = _SYNTHETIC_FIELD_CONSTRUCTORS[func]
                    synth_function(grid, name, at, **name_dict)
                elif func is "read_esri_ascii":
                    read_esri_ascii(grid=grid, **name_dict)
                elif func is "read_netcdf":
                    read_netcdf(grid=grid, **name_dict)
                else:
                    msg = "Bad function supplied to construct field"
                    raise ValueError(msg)

    if len(fields_dict) is not 0:
        msg = "Bad group location supplied to construct fields."
        raise ValueError(msg)

    # part three, set boundary conditions
    bc_list = dict_like.pop("boundary_conditions")
    while len(bc_list) > 0:
        bc_function = bc_list.pop(0)
        if bc_function in grid.__dict__:
            pass

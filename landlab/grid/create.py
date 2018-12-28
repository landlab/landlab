#! /usr/bin/env python
"""Create landlab model grids."""
from warnings import warn

from landlab.core import model_parameter_dictionary as mpd
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
    """Create grid, initialize fields, and set boundary conditions.

    Examples
    --------
    >>> import numpy as np
    >>> np.random.seed(42)
    >>> from landlab import create_grid
    >>> p = {'grid': {'RasterModelGrid': {'shape': (4,5),
    ...                                   'xy_spacing': (3, 4)}
    ...               },
    ...      'fields': {'at_node': {'spam': {'plane': {'point': (1, 1, 1),
    ...                                                'normal': (-2, -1, 1)},
    ...                                      'random': {'distribution': 'uniform',
    ...                                                 'low': 1,
    ...                                                 'high': 4}
    ...                                      },
    ...                             },
    ...                 'at_link': {'eggs' : {'constant': {'where': 'ACTIVE_LINK',
    ...                                                    'constant': 12},
    ...                                       },
    ...                             },
    ...                 },
    ...      'boundary_conditions': [
    ...                     {'set_closed_boundaries_at_grid_edges':
    ...                                 {'right_is_closed': True,
    ...                                  'top_is_closed': True,
    ...                                  'left_is_closed': True,
    ...                                  'bottom_is_closed': True
    ...                                  }
    ...                      }]
    ...      }
    >>> mg = create_grid(p)
    >>> mg.number_of_nodes
    20
    >>> "spam" in mg.at_node
    True
    >>> "eggs" in mg.at_link
    True

    """
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

    grid_params = grid_dict.pop(grid_type)
    grid = grid_class.from_dict(grid_params)

    # part two, create fields
    fields_dict = dict_like.pop("fields", {})

    # for each grid element:
    for at_group in fields_dict:
        at = at_group[3:]
        if at not in grid.groups:
            msg = ("")
            raise ValueError(msg)

        at_dict = fields_dict[at_group]
        # for field at grid element
        for name in at_dict:
            name_dict = at_dict[name]

            # for each function, add values.
            for func in name_dict:
                func_dict = name_dict[func]
                if func in _SYNTHETIC_FIELD_CONSTRUCTORS:
                    synth_function = _SYNTHETIC_FIELD_CONSTRUCTORS[func]
                    synth_function(grid, name, at=at, **func_dict)
                elif func is "read_esri_ascii":
                    read_esri_ascii(grid=grid, **func_dict)
                elif func is "read_netcdf":
                    read_netcdf(grid=grid, **func_dict)
                else:
                    msg = "Bad function supplied to construct field"

    # part three, set boundary conditions
    bc_list = dict_like.pop("boundary_conditions")
    for bc_function_dict in bc_list:

        if len(bc_function_dict) != 1:
            msg = ""
            raise ValueError(msg)

        for bc_function in bc_function_dict:
            bc_params = bc_function_dict[bc_function]

            # these two will work for all types.
            if bc_function is "set_nodata_nodes_to_closed":
                grid.set_nodata_nodes_to_closed(**bc_params)
            elif bc_function is "set_nodata_nodes_to_fixed_gradient":
                grid.set_nodata_nodes_to_fixed_gradient(**bc_params)

            # this set only works for raster
            elif (
                bc_function is "set_closed_boundaries_at_grid_edges"
                and grid_type is "RasterModelGrid"
            ):
                grid.set_closed_boundaries_at_grid_edges(**bc_params)
            elif (
                bc_function is "set_fixed_link_boundaries_at_grid_edges"
                and grid_type is "RasterModelGrid"
            ):
                grid.set_fixed_link_boundaries_at_grid_edges(**bc_params)
            elif (
                bc_function is "set_fixed_value_boundaries_at_grid_edges"
                and grid_type is "RasterModelGrid"
            ):
                grid.set_fixed_value_boundaries_at_grid_edges(**bc_params)
            elif (
                bc_function is "set_watershed_boundary_condition_outlet_coords"
                and grid_type is "RasterModelGrid"
            ):
                grid.set_watershed_boundary_condition_outlet_coords(**bc_params)
            elif (
                bc_function is "set_open_nodes_disconnected_from_watershed_to_closed"
                and grid_type is "RasterModelGrid"
            ):
                grid.set_open_nodes_disconnected_from_watershed_to_closed(**bc_params)
            elif (
                bc_function is "set_status_at_node_on_edges"
                and grid_type is "RasterModelGrid"
            ):
                grid.set_status_at_node_on_edges(**bc_params)

            # these work for hex and raster only.
            elif (
                bc_function is "set_open_nodes_disconnected_from_watershed_to_closed"
                and grid_type in ("RasterModelGrid", "HexModelGrid")
            ):
                grid.set_watershed_boundary_condition(**bc_params)
            elif (
                bc_function is "set_open_nodes_disconnected_from_watershed_to_closed"
                and grid_type in ("RasterModelGrid", "HexModelGrid")
            ):
                grid.set_watershed_boundary_condition(**bc_params)

    return grid

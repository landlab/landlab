#! /usr/bin/env python
"""Create landlab model grids."""


from landlab.core import model_parameter_dictionary as mpd
from .raster import from_dict as raster_from_dict
from .hex import from_dict as hex_from_dict


class Error(Exception):

    """Base class for exceptions from this module."""

    pass


class BadGridTypeError(Error):

    """Raise this error for a bad grid type."""

    def __init__(self, grid_type):
        self._type = str(grid_type)

    def __str__(self):
        return self._type


_GRID_READERS = {
    'raster': raster_from_dict,
    'hex': hex_from_dict,
}


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
    >>> grid = create_and_initialize_grid(test_file)
    >>> grid.number_of_nodes
    20
    """
    if isinstance(input_source, dict):
        param_dict = input_source
    else:
        param_dict = mpd.ModelParameterDictionary(from_file=input_source)

    grid_type = param_dict['GRID_TYPE']

    grid_type.strip().lower()

    # Read parameters appropriate to that type, create it, and initialize it
    try:
        grid_reader = _GRID_READERS[grid_type]
    except KeyError:
        raise BadGridTypeError(grid_type)

    # Return the created and initialized grid
    return grid_reader(param_dict)

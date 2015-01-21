#! /usr/bin/env python

import types

from landlab.core import model_parameter_dictionary as mpd
from .raster import from_dict as raster_from_dict
from .hex import from_dict as hex_from_dict


class Error(Exception):
    pass


class BadGridTypeError(Error):
    def __init__(self, grid_type):
        self._type = str(grid_type)

    def __str__(self):
        return self._type


_GRID_READERS = {
    'raster': raster_from_dict,
    'hex': hex_from_dict,
}


def create_and_initialize_grid(input_source):
    """
    Creates, initializes, and returns a new grid object using parameters 
    specified in *input_source*. *input_source* is either a
    ModelParameterDictionary instance (or, really, just dict-like) or a
    named input file.
    
    Example:
        
    >>> from StringIO import StringIO
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
    >>> mg = create_and_initialize_grid(test_file)
    >>> mg.number_of_nodes
    20
        
    """    
    # Handle input source. 
    #In this code block, we do the following:
    #   - handle the case in which caller provides neither a parameter
    #     dictionary nor an input file name.
    #   - if we're given an input file name, create a parameter dictionary
    #     object that reads the specified file name
    if isinstance(input_source, dict):
        param_dict = input_source
    else:
        param_dict = mpd.ModelParameterDictionary(from_file=input_source)
        
    # Find out what type of grid the user wants
    #
    # Dev note: could handle defaults like:
    #    param_dict.get('GRID_TYPE','raster')
    # so if no GRID_TYPE is specified you get the second arg as default. If no
    # second arg, then exception
    try:
        grid_type = param_dict['GRID_TYPE']
    except KeyError:
        raise

    grid_type.strip().lower()   # make LC w/o leading/trailing spaces
    
    # Read parameters appropriate to that type, create it, and initialize it
    try:
        grid_reader = _GRID_READERS[grid_type]
    except KeyError:
        raise BadGridTypeError(grid_type)
    else:
        mg = grid_reader(param_dict)
    
    # Return the created and initialized grid
    return mg

#! /usr/bin/env python
"""
Read data from an ESRI ASCII file into a RasterModelGrid
"""

import os
import types
import re

import numpy as np

import landlab
from landlab import RasterModelGrid


_VALID_HEADER_KEYS = [
    'ncols', 'nrows', 'xllcorner', 'xllcenter', 'yllcorner',
    'yllcenter', 'cellsize', 'nodata_value',
]
_HEADER_REGEX_PATTERN = re.compile('\s*(?P<key>\w+)\s+(?P<value>[\w.+-]+)')
_HEADER_VALUE_TESTS = {
    'nrows': (int, lambda x: x > 0),
    'ncols': (int, lambda x: x > 0),
    'cellsize': (float, lambda x: x > 0),
    'xllcorner': (float, lambda x: True),
    'xllcenter': (float, lambda x: True),
    'yllcorner': (float, lambda x: True),
    'yllcenter': (float, lambda x: True),
    'nodata_value': (float, lambda x: True),
}


class Error(Exception):
    """
    Base class for errors in this module.
    """
    pass


class BadHeaderLineError(Error):
    """
    Raise this error if a bad header is line encounter in a ESRI ASCII
    file.
    """
    def __init__(self, line):
        self._line = line

    def __str__(self):
        return self._line


class MissingRequiredKeyError(Error):
    """
    Raise this error when a header is missing a required key.
    """
    def __init__(self, key):
        self._key = key

    def __str__(self):
        return self._key


class KeyTypeError(Error):
    """
    Raise this error when a header's key value is of the wrong type.
    """
    def __init__(self, key, expected_type):
        self._key = key
        self._type = str(expected_type)

    def __str__(self):
        return 'Unable to convert %s to %s' % (self._key, self._type)


class KeyValueError(Error):
    """
    Raise this error when a header's key value has a bad value
    """
    def __init__(self, key, message):
        self._key = key
        self._msg = message

    def __str__(self):
        return '%s: %s' % (self._key, self._msg)


class DataSizeError(Error):
    """
    Raise this error if the size of the data does not match that given in
    the header.
    """
    def __init__(self, size, expected_size):
        self._actual = size
        self._expected = expected_size

    def __str__(self):
        return '%s != %s' % (self._actual, self._expected)


def _parse_header_key_value(line):
    match = _HEADER_REGEX_PATTERN.match(line)
    if match is None:
        raise BadHeaderLineError(line)

    (key, value) = (match.group('key').lower(), match.group('value'))
    try:
        assert(key in _VALID_HEADER_KEYS)
    except AssertionError:
        raise BadHeaderLineError(line)
    else:
        return (key, value)


def _header_lines(asc_file):
    pos = asc_file.tell()
    line = asc_file.readline()
    while len(line) > 0:
        if len(line.strip()) > 0:
            try:
                (key, value) = _parse_header_key_value(line)
            except BadHeaderLineError:
                asc_file.seek(pos, 0)
                break
            else:
                yield (key, value)
        pos = asc_file.tell()
        line = asc_file.readline()


def _header_is_valid(header):
    header_keys = set(header)
    required_keys = set(['ncols', 'nrows', 'cellsize'])

    try:
        assert(required_keys.issubset(header_keys))
    except AssertionError:
        raise MissingRequiredKeyError(', '.join(required_keys - header_keys))

    for keys in [('xllcenter', 'xllcorner'), ('yllcenter', 'yllcorner')]:
        try:
            assert(len(set(keys) & header_keys) == 1)
        except AssertionError:
            raise MissingRequiredKeyError('|'.join(keys))

    for (key, requires) in _HEADER_VALUE_TESTS.items():
        try:
            header[key] = requires[0](header[key])
            assert(requires[1](header[key]))
        except ValueError:
            raise KeyTypeError(key, float)
        except AssertionError:
            raise KeyValueError(key, 'Bad value')
        except KeyError:
            pass

    return True

def read_asc_header(asc_file):
    """
    Read header information from an ESRI ASCII raster file.

    The header contains the following variables,
        - *ncols*: Number of cell columns
        - *nrows*: Number of cell rows
        - *xllcenter* or *xllcorner*: X (column) coordinate of lower-left
            coordinate of grid (by center or lower-left corner of the cell)
        - *yllcenter*, *yllcorner*: Y (row) coordinate of lower-left
            coordinate of grid (by center or lower-left corner of the cell)
        - *cellsize*: Grid spacing between rows and columns
        - *nodata_value*: No-data value (optional)
    """
    header = dict()
    for (key, value) in _header_lines(asc_file):
        header[key] = value

    try:
        _header_is_valid(header)
    except Error:
        raise

    return header


def _read_asc_data(asc_file, header={}):
    """
    Read gridded data from an ESRI ASCII data file.

    :asc_file: File-like object of the data file pointing to the start of the
               data.

    .. note::
        Row 1 of the data is at the top of the raster, row 2 is just under
        row 1, and so on.
    """
    #return np.fromtxt(asc_file, skiprows=len(header))
    try:
        return np.loadtxt(asc_file)
    except ValueError as error:
        print header
        print error
    #return np.genfromtxt(asc_file)


def read_esri_ascii(asc_file, reshape=False):
    """Read :py:class:`~landlab.RasterModelGrid` from an ESRI ASCII file.

    Read data from *asc_file*, an ESRI_ ASCII file, into a
    :py:class:`~landlab.RasterModelGrid`.  *asc_file* is either the name of
    the data file or is a file-like object.

    The grid and data read from the file are returned as a tuple
    (*grid*, *data*) where *grid* is an instance of
    :py:class:`~landlab.RasterModelGrid` and *data* is a numpy
    array of doubles with that has been reshaped to have the number of rows
    and columns given in the header.

    .. _ESRI: http://resources.esri.com/help/9.3/arcgisengine/java/GP_ToolRef/spatial_analyst_tools/esri_ascii_raster_format.htm

    Parameters
    ----------
    asc_file : str of file-like
        Data file to read.
    reshape : boolean, optional
        Reshape the returned array, otherwise return a flattened array.

    Returns
    -------
    (grid, data) : tuple
        A newly-created RasterModel grid and the associated node data.
    """
    if isinstance(asc_file, types.StringTypes):
        file_name = asc_file
        with open(file_name, 'r') as asc_file:
            header = read_asc_header(asc_file)
            data = _read_asc_data(asc_file, header=header)
    else:
        header = read_asc_header(asc_file)
        data = _read_asc_data(asc_file, header=header)

    shape = (header['nrows'], header['ncols'])
    spacing = (header['cellsize'], header['cellsize'])
    origin = (header['xllcorner'], header['yllcorner'])

    try:
        assert(data.size == shape[0] * shape[1])
    except AssertionError:
        raise DataSizeError(shape[0] * shape[1], data.size)

    data.shape = shape
    data = np.flipud(data)
    if not reshape:
        data = data.flatten()

    grid = RasterModelGrid(num_rows=shape[0], num_cols=shape[1],
                           dx=spacing[0])

    return (grid, data)


def write_esri_ascii(path, fields, names=None, clobber=False):
    """Write landlab fields to ESRI ASCII.

    Write the data and grid information for *fields* to *path* in the ESRI
    ASCII format.

    Parameters
    ----------
    path : str
        Path to output file.
    fields : field-like
        Landlab field object that holds a grid and associated values.
    names : iterable of str, optional
        Names of the fields to include in the output file. If not provided,
        write all fields.
    clobber : boolean
        If *path* exists, clobber the existing file, otherwise raise an
        exception.

    Examples
    --------
    >>> import numpy as np
    >>> from landlab.testing.tools import cdtemp
    >>> from landlab import RasterModelGrid
    >>> from landlab.io.esri_ascii import write_esri_ascii

    >>> grid = RasterModelGrid(4, 5, dx=2.)
    >>> _ = grid.add_field('node', 'air__temperature', np.arange(20.))
    >>> with cdtemp() as _:
    ...     files = write_esri_ascii('test.asc', grid)
    >>> files
    ['test.asc']

    >>> _ = grid.add_field('node', 'land_surface__elevation', np.arange(20.))
    >>> with cdtemp() as _:
    ...     files = write_esri_ascii('test.asc', grid)
    >>> files.sort()
    >>> files
    ['test_air__temperature.asc', 'test_land_surface__elevation.asc']
    """
    if os.path.exists(path) and not clobber:
        raise ValueError('file exists')

    if isinstance(names, types.StringTypes):
        names = [names]

    names = names or fields.at_node.keys()

    if len(names) == 1:
        paths = [path]
    elif len(names) > 1:
        (base, ext) = os.path.splitext(path)
        paths = [base + '_' + name + ext for name in names]
    else:
        raise ValueError('no node fields to write')

    bad_names = set(names) - set(fields.at_node.keys())
    if len(bad_names) > 0:
        raise ValueError('unknown field name(s): %s' % ','.join(bad_names))

    header = {
        'ncols': fields.number_of_node_columns,
        'nrows': fields.number_of_node_rows,
        'xllcorner': fields.node_x[0],
        'yllcorner': fields.node_y[0],
        'cellsize': fields.dx,
    }

    for path, name in zip(paths, names):
        with open(path, 'w') as fp:
            for key, val in header.items():
                fp.write('%s %s%s' % (key, str(val), os.linesep))
            data = fields.at_node[name].reshape(header['nrows'],
                                                header['ncols'])
            np.savetxt(fp, np.flipud(data))

    return paths


def write_esri_ascii_from_raster_model_grid_data(grid, data_name=None,
                                                 data_array=None, xllc=0.,
                                                 yllc=0., nodata_value=-9999,
                                                 filename='landlab_esri_ascii_grid.asc'):
    """
    Creates a new ESRI-Ascii-grid format file(s) containing data attached to a
    RasterModelGrid.

    Parameters
    ----------

    grid : RasterModelGrid
        Grid on which data are organized
    data_name : string
        Name of data item belonging to grid (must give either this or data_array)
    data_array : Numpy 1D array with data corresponding to nodes
        Data to be output (must give either this or data_name)
    xllc, yllc : float or int
        x and y coordinates of grid's lower-left corner
    nodata_value : int or float
        Code used to represent a grid cell containing invalid/no data
    filename : string
        Name of file to create


    Returns
    -------
    None
    """

    # Make sure we've been passed a RasterModelGrid
    assert (type(grid) is landlab.grid.raster.RasterModelGrid), \
           '<grid> must be a Landlab RasterModelGrid'

    # Make sure user has specified either a data-item name or an array
    assert ((data_name is not None) or (data_array is not None)), \
            'Must pass a non-None value for either data_name or data_array'

    # If data_name specified, find data with that name at nodes. If not, then
    # use data_array. If that's not specified either, throw an exception.
    if data_name is not None:
        data = grid.at_node[data_name][::-1,:]
    else:
        data = data_array[::-1,:]

    # Open a new text file
    output_file = open(filename, 'w')

    # Write the header
    output_file.write('ncols        %d\n' % (grid.number_of_node_columns))
    output_file.write('nrows        %d\n' % (grid.number_of_node_rows))
    output_file.write('xllcorner    %.3f\n' % (xllc))
    output_file.write('yllcorner    %.3f\n' % (yllc))
    output_file.write('cellsize     %.3f\n' % (grid.dx))
    output_file.write('NODATA_value %.0f\n' % (nodata_value))

    # Write the data, row by row
    n = 0
    for r in range(grid.number_of_node_rows):
        for c in range(grid.number_of_node_columns):
            output_file.write('%.3f ' % (data[n]))
            n+=1
        output_file.write('\n')

    # Close the file
    output_file.close()

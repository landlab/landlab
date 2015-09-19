#! /usr/bin/env python
"""Read data from an ESRI ASCII file into a RasterModelGrid."""

import os
import re
import six

import numpy as np

_VALID_HEADER_KEYS = [
    'ncols', 'nrows', 'xllcorner', 'xllcenter', 'yllcorner',
    'yllcenter', 'cellsize', 'nodata_value',
]
_HEADER_KEY_REGEX_PATTERN = re.compile(r'\s*(?P<key>[a-zA-z]\w+)')
_HEADER_REGEX_PATTERN = re.compile(
    r'\s*(?P<key>[a-zA-Z]\w+)\s+(?P<value>[\w.+-]+)')
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

    """Base class for errors in this module."""

    pass


class BadHeaderLineError(Error):

    """Raise this error for a bad header is line."""

    def __init__(self, line):
        self._line = line

    def __str__(self):
        return self._line


class MissingRequiredKeyError(Error):

    """Raise this error when a header is missing a required key."""

    def __init__(self, key):
        self._key = key

    def __str__(self):
        return self._key


class KeyTypeError(Error):

    """Raise this error when a header's key value is of the wrong type."""

    def __init__(self, key, expected_type):
        self._key = key
        self._type = str(expected_type)

    def __str__(self):
        return 'Unable to convert %s to %s' % (self._key, self._type)


class KeyValueError(Error):

    """Raise this error when a header's key value has a bad value."""

    def __init__(self, key, message):
        self._key = key
        self._msg = message

    def __str__(self):
        return '%s: %s' % (self._key, self._msg)


class DataSizeError(Error):

    """Raise this error if the size of data does not match the header."""

    def __init__(self, size, expected_size):
        self._actual = size
        self._expected = expected_size

    def __str__(self):
        return '%s != %s' % (self._actual, self._expected)


def _parse_header_key_value(line):
    """Parse a header line into a key-value pair.

    Parameters
    ----------
    line : str
        Header line.

    Returns
    -------
    (str, str)
        Header key-value pair

    Raises
    ------
    BadHeaderLineError
        The is something wrong with the header line.
    """
    match = _HEADER_KEY_REGEX_PATTERN.match(line)
    if match is None:
        return None
        # raise BadHeaderLineError(line)

    match = _HEADER_REGEX_PATTERN.match(line)
    if match is None:
        raise BadHeaderLineError(line)

    (key, value) = (match.group('key').lower(), match.group('value'))

    if key in _VALID_HEADER_KEYS:
        return (key, value)
    else:
        raise BadHeaderLineError(line)


def _header_lines(asc_file):
    """Iterate over header lines for a ESRI ASCII file.

    Parameters
    ----------
    asc_file : file_like
        File-like object for an ESRI ASCII file.

    Yields
    ------
    str
        Header line.
    """
    pos = asc_file.tell()
    line = asc_file.readline()
    while len(line) > 0:
        if len(line.strip()) > 0:
            item = _parse_header_key_value(line)
            if item:
                yield item
            else:
                asc_file.seek(pos, 0)
                break
        pos = asc_file.tell()
        line = asc_file.readline()


def _header_is_valid(header):
    """Check if the ESRI ASCII header is valid.

    Parameters
    ----------
    header : dict
        Header as key-values pairs.

    Raises
    ------
    MissingRequiredKeyError
        The header is missing a required key.
    KeyTypeError
        The header has the key but its values is of the wrong type.
    """
    header_keys = set(header)
    required_keys = set(['ncols', 'nrows', 'cellsize'])

    if not required_keys.issubset(header_keys):
        raise MissingRequiredKeyError(', '.join(required_keys - header_keys))

    for keys in [('xllcenter', 'xllcorner'), ('yllcenter', 'yllcorner')]:
        if len(set(keys) & header_keys) != 1:
            raise MissingRequiredKeyError('|'.join(keys))

    for (key, requires) in _HEADER_VALUE_TESTS.items():
        to_type, is_valid = requires

        if key not in header:
            continue

        try:
            header[key] = to_type(header[key])
        except ValueError:
            raise KeyTypeError(key, to_type)

        if not is_valid(header[key]):
            raise KeyValueError(key, 'Bad value')

    return True


def read_asc_header(asc_file):
    """Read header information from an ESRI ASCII raster file.

    The header contains the following variables,
        - *ncols*: Number of cell columns
        - *nrows*: Number of cell rows
        - *xllcenter* or *xllcorner*: X (column) coordinate of lower-left
            coordinate of grid (by center or lower-left corner of the cell)
        - *yllcenter*, *yllcorner*: Y (row) coordinate of lower-left
            coordinate of grid (by center or lower-left corner of the cell)
        - *cellsize*: Grid spacing between rows and columns
        - *nodata_value*: No-data value (optional)

    Parameters
    ----------
    asc_file : file_like
        File-like object from which to read header.

    Returns
    -------
    dict
        Header as key-value pairs.

    Raises
    ------
    MissingRequiredKeyError
        The header is missing a required key.
    KeyTypeError
        The header has the key but its values is of the wrong type.

    Examples
    --------
    >>> from six import StringIO
    >>> from landlab.io.esri_ascii import read_asc_header
    >>> contents = StringIO('''
    ...     nrows 100
    ...     ncols 200
    ...     cellsize 1.5
    ...     xllcenter 0.5
    ...     yllcenter -0.5
    ... ''')
    >>> hdr = read_asc_header(contents)
    >>> hdr['nrows'], hdr['ncols']
    (100, 200)
    >>> hdr['cellsize']
    1.5
    >>> hdr['xllcenter'], hdr['yllcenter']
    (0.5, -0.5)

    ``MissingRequiredKey`` is raised if the header does not contain all of the
    necessary keys.

    >>> contents = StringIO('''
    ...     ncols 200
    ...     cellsize 1.5
    ...     xllcenter 0.5
    ...     yllcenter -0.5
    ... ''')
    >>> read_asc_header(contents)
    Traceback (most recent call last):
        ...
    MissingRequiredKeyError: nrows

    ``KeyTypeError`` is raises if a value is of the wrong type. For instance,
    ``nrows`` and ``ncols`` must be ``int``.

    >>> contents = StringIO('''
    ...     nrows 100.5
    ...     ncols 200
    ...     cellsize 1.5
    ...     xllcenter 0.5
    ...     yllcenter -0.5
    ... ''')
    >>> read_asc_header(contents)
    Traceback (most recent call last):
        ...
    KeyTypeError: Unable to convert nrows to <type 'int'>
    """
    header = dict()
    for (key, value) in _header_lines(asc_file):
        header[key] = value

    _header_is_valid(header)

    return header


def _read_asc_data(asc_file):
    """Read gridded data from an ESRI ASCII data file.

    Parameters
    ----------
    asc_file : file-like
        File-like object of the data file pointing to the start of the data.

    .. note::
        First row of the data is at the top of the raster grid, the second
        row is the second from the top, and so on.
    """
    return np.loadtxt(asc_file)


def read_esri_ascii(asc_file, grid=None, reshape=False, name=None):
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
    name : str, optional
        Add data to the grid as a named field.
    grid : *grid* , optional
        Adds data to an existing *grid* instead of creating a new one.

    Returns
    -------
    (grid, data) : tuple
        A newly-created RasterModel grid and the associated node data.
    """
    from ..grid import RasterModelGrid

    if isinstance(asc_file, six.string_types):
        file_name = asc_file
        with open(file_name, 'r') as asc_file:
            header = read_asc_header(asc_file)
            data = _read_asc_data(asc_file)
    else:
        header = read_asc_header(asc_file)
        data = _read_asc_data(asc_file)

    shape = (header['nrows'], header['ncols'])
    spacing = (header['cellsize'], header['cellsize'])
    origin = (header['xllcorner'], header['yllcorner'])

    if data.size != shape[0] * shape[1]:
        raise DataSizeError(shape[0] * shape[1], data.size)

    data.shape = shape
    data = np.flipud(data)
    if not reshape:
        data = data.flatten()

    if grid is None:
        grid = RasterModelGrid(num_rows=shape[0], num_cols=shape[1],
                               dx=spacing[0])
    if name:
        grid.add_field('node', name, data)

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

    if isinstance(names, six.string_types):
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
        header_lines = ['%s %s' % (key, str(val))
                        for key, val in list(header.items())]
        data = fields.at_node[name].reshape(header['nrows'], header['ncols'])
        np.savetxt(path, np.flipud(data), header=os.linesep.join(header_lines),
                   comments='')

    return paths

#! /usr/bin/env python
"""Read/write data from an ESRI ASCII file into a RasterModelGrid.

ESRI ASCII functions
++++++++++++++++++++

.. autosummary::

    ~landlab.io.esri_ascii.read_asc_header
    ~landlab.io.esri_ascii.read_esri_ascii
    ~landlab.io.esri_ascii.write_esri_ascii
"""

import os
import pathlib
import re

import numpy as np

from landlab.utils import add_halo

_VALID_HEADER_KEYS = [
    "ncols",
    "nrows",
    "xllcorner",
    "xllcenter",
    "yllcorner",
    "yllcenter",
    "cellsize",
    "nodata_value",
]
_HEADER_KEY_REGEX_PATTERN = re.compile(r"\s*(?P<key>[a-zA-z]\w+)")
_HEADER_REGEX_PATTERN = re.compile(r"\s*(?P<key>[a-zA-Z]\w+)\s+(?P<value>[\w.+-]+)")
_HEADER_VALUE_TESTS = {
    "nrows": (int, lambda x: x > 0),
    "ncols": (int, lambda x: x > 0),
    "cellsize": (float, lambda x: x > 0),
    "xllcorner": (float, lambda x: True),
    "xllcenter": (float, lambda x: True),
    "yllcorner": (float, lambda x: True),
    "yllcenter": (float, lambda x: True),
    "nodata_value": (float, lambda x: True),
}


class Error(Exception):
    """Base class for errors in this module."""

    pass


class BadHeaderLineError(Error):
    """Raise this error for a bad header is line."""

    def __init__(self, line):
        self._line = line

    def __str__(self):
        return self._line  # this line not yet tested


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
        return f"Unable to convert {self._key} to {self._type}"


class KeyValueError(Error):
    """Raise this error when a header's key value has a bad value."""

    def __init__(self, key, message):
        self._key = key
        self._msg = message

    def __str__(self):
        return f"{self._key}: {self._msg}"  # this line not yet tested


class DataSizeError(Error):
    """Raise this error if the size of data does not match the header."""

    def __init__(self, size, expected_size):
        self._actual = size
        self._expected = expected_size

    def __str__(self):
        return "{} != {}".format(
            self._actual, self._expected
        )  # this line not yet tested


class MismatchGridDataSizeError(Error):
    """Raise this error if the data size does not match the grid size."""

    def __init__(self, size, expected_size):
        self._actual = size
        self._expected = expected_size

    def __str__(self):
        return "(data size) {} != {} (grid size)".format(
            self._actual,
            self._expected,
        )  # this line not yet tested


class MismatchGridXYSpacing(Error):
    """Raise this error if the file cell size does not match the grid dx."""

    def __init__(self, dx, expected_dx):
        self._actual = dx
        self._expected = expected_dx

    def __str__(self):
        return "(data dx) {} != {} (grid dx)".format(
            self._actual,
            self._expected,
        )  # this line not yet tested


class MismatchGridXYLowerLeft(Error):
    """Raise this error if the file lower left does not match the grid."""

    def __init__(self, llc, expected_llc):
        self._actual = llc
        self._expected = expected_llc

    def __str__(self):
        return "(data lower-left) {} != {} (grid lower-left)".format(
            self._actual,
            self._expected,
        )  # this line not yet tested


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

    (key, value) = (match.group("key").lower(), match.group("value"))

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
    required_keys = {"ncols", "nrows", "cellsize"}

    if not required_keys.issubset(header_keys):
        raise MissingRequiredKeyError(", ".join(required_keys - header_keys))

    for keys in [("xllcenter", "xllcorner"), ("yllcenter", "yllcorner")]:
        if len(set(keys) & header_keys) != 1:
            raise MissingRequiredKeyError("|".join(keys))

    for key, requires in _HEADER_VALUE_TESTS.items():
        to_type, is_valid = requires

        if key not in header:
            continue

        try:
            header[key] = to_type(header[key])
        except ValueError as exc:
            raise KeyTypeError(key, to_type) from exc

        if not is_valid(header[key]):
            raise KeyValueError(key, "Bad value")

    return True


def read_asc_header(asc_file):
    """Read header information from an ESRI ASCII raster file.

    The header contains the following variables,

    * ``ncols``: Number of cell columns
    * ``nrows``: Number of cell rows
    * ``xllcenter`` or ``xllcorner``: X (column) coordinate of lower-left
        coordinate of grid (by center or lower-left corner of the cell)
    * ``yllcenter``, ``yllcorner``: Y (row) coordinate of lower-left
        coordinate of grid (by center or lower-left corner of the cell)
    * ``cellsize``: Grid spacing between rows and columns
    * ``nodata_value``: No-data value (optional)

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
    :class:`~landlab.io.esri_ascii.MissingRequiredKeyError`
        The header is missing a required key.
    :class:`~landlab.io.esri_ascii.KeyTypeError`
        The header has the key but its values is of the wrong type.

    Examples
    --------
    >>> from io import StringIO
    >>> from landlab.io.esri_ascii import read_asc_header

    >>> contents = '''
    ... nrows 100
    ... ncols 200
    ... cellsize 1.5
    ... xllcenter 0.5
    ... yllcenter -0.5
    ... '''

    >>> hdr = read_asc_header(StringIO(contents))
    >>> hdr["nrows"], hdr["ncols"]
    (100, 200)
    >>> hdr["cellsize"]
    1.5
    >>> hdr["xllcenter"], hdr["yllcenter"]
    (0.5, -0.5)

    :class:`~landlab.io.esri_ascii.MissingRequiredKeyError` is raised if the
    header does not contain all of the necessary keys.

    >>> contents = '''
    ... ncols 200
    ... cellsize 1.5
    ... xllcenter 0.5
    ... yllcenter -0.5
    ... '''
    >>> read_asc_header(StringIO(contents))
    Traceback (most recent call last):
    MissingRequiredKeyError: nrows

    :class:`~landlab.io.esri_ascii.KeyTypeError` is raised if a value is of
    the wrong type. For instance, *nrows* and *ncols* must be ``int``.

    >>> contents = '''
    ... nrows 100.5
    ... ncols 200
    ... cellsize 1.5
    ... xllcenter 0.5
    ... yllcenter -0.5
    ... '''
    >>> read_asc_header(StringIO(contents))
    Traceback (most recent call last):
    KeyTypeError: Unable to convert nrows to <type 'int'>
    """
    header = {}
    for key, value in _header_lines(asc_file):
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


def read_esri_ascii(asc_file, grid=None, reshape=False, name=None, halo=0):
    """Read :py:class:`~.RasterModelGrid` from an ESRI ASCII file.

    Read data from *asc_file*, an `ESRI ASCII file`_, into a
    :class:`~.RasterModelGrid`.  *asc_file* is either the name of
    the data file or is a file-like object.

    The grid and data read from the file are returned as a tuple
    (*grid*, *data*) where *grid* is an instance of
    :py:class:`~.RasterModelGrid` and *data* is a numpy
    array of doubles with that has been reshaped to have the number of rows
    and columns given in the header.

    .. _ESRI ASCII file: http://resources.esri.com/help/9.3/arcgisengine/java/GP_ToolRef/spatial_analyst_tools/esri_ascii_raster_format.htm

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
    halo : integer, optional
        Adds outer border of depth halo to the *grid*.

    Returns
    -------
    (grid, data) : tuple
        A newly-created :class:`~.RasterModelGrid` and the associated node data.

    Raises
    ------
    :class:`~landlab.io.esri_ascii.DataSizeError`
        Data are not the same size as indicated by the header file.
    :class:`~landlab.io.esri_ascii.MismatchGridDataSizeError`
        If a grid is passed, and the size of the grid does not agree with the
        size of the data.
    :class:`~landlab.io.esri_ascii.MismatchGridXYSpacing`
        If a grid is passed, and the *cellsize* listed in the heading does not
        match the node spacing of the grid.
    :class:`~landlab.io.esri_ascii.MismatchGridXYLowerLeft`
        If a grid is passed and the *xllcorner* and *yllcorner* do not match that
        of the grid.

    Examples
    --------

    >>> from landlab.io import read_esri_ascii
    >>> from io import StringIO

    >>> contents = '''
    ... ncols         3
    ... nrows         4
    ... xllcorner     1.
    ... yllcorner     2.
    ... cellsize      10.
    ... NODATA_value  -1
    ... 0. 1. 2.
    ... 3. 4. 5.
    ... 6. 7. 8.
    ... 9. 10. 11.
    ... '''

    >>> (grid, data) = read_esri_ascii(StringIO(contents))

    The returned grid is a :class:`~.RasterModelGrid` with 4 rows and 3 columns.

    >>> grid
    RasterModelGrid((4, 3), xy_spacing=(10.0, 10.0), xy_of_lower_left=(1.0, 2.0))

    Note that the first row of values is the bottom-most of the data file.

    >>> data.reshape(grid.shape)
    array([[  9.,  10.,  11.],
           [  6.,   7.,   8.],
           [  3.,   4.,   5.],
           [  0.,   1.,   2.]])

    >>> (grid, data) = read_esri_ascii(StringIO(contents), halo=1)

    Because of the halo, the returned grid now has two more rows and columns than before.

    >>> grid
    RasterModelGrid((6, 5), xy_spacing=(10.0, 10.0), xy_of_lower_left=(-9.0, -8.0))
    >>> data.reshape(grid.shape)
    array([[-1.,  -1.,  -1.,  -1.,  -1.],
           [-1.,   9.,  10.,  11.,  -1.],
           [-1.,   6.,   7.,   8.,  -1.],
           [-1.,   3.,   4.,   5.,  -1.],
           [-1.,   0.,   1.,   2.,  -1.],
           [-1.,  -1.,  -1.,  -1.,  -1.]])
    """  # noqa: B950
    from ..grid import RasterModelGrid

    if halo < 0:
        raise ValueError("negative halo")

    # if the asc_file is provided as a string, open it and pass the pointer to
    # _read_asc_header, and _read_asc_data
    if isinstance(asc_file, (str, pathlib.Path)):
        with open(asc_file) as f:
            header = read_asc_header(f)
            data = _read_asc_data(f)

    # otherwise, pass asc_file directly.
    else:
        header = read_asc_header(asc_file)
        data = _read_asc_data(asc_file)

    shape = (header["nrows"] + 2 * halo, header["ncols"] + 2 * halo)
    nodata_value = header.get("nodata_value", -9999.0)
    if data.size != (shape[0] - 2 * halo) * (shape[1] - 2 * halo):
        raise DataSizeError(shape[0] * shape[1], data.size)

    xy_spacing = (header["cellsize"], header["cellsize"])
    xy_of_lower_left = (
        header["xllcorner"] - halo * header["cellsize"],
        header["yllcorner"] - halo * header["cellsize"],
    )

    data = np.flipud(data)

    if halo > 0:
        data = add_halo(
            data.reshape(header["nrows"], header["ncols"]),
            halo=halo,
            halo_value=nodata_value,
        ).reshape((-1,))

    if not reshape:
        data = data.flatten()

    if grid is not None:
        if (grid.number_of_node_rows != shape[0]) or (
            grid.number_of_node_columns != shape[1]
        ):
            raise MismatchGridDataSizeError(
                shape[0] * shape[1],
                grid.number_of_node_rows * grid.number_of_node_columns,
            )
        if (grid.dx, grid.dy) != xy_spacing:
            raise MismatchGridXYSpacing((grid.dx, grid.dy), xy_spacing)

        if grid.xy_of_lower_left != xy_of_lower_left:
            raise MismatchGridXYLowerLeft(grid.xy_of_lower_left, xy_of_lower_left)

    if grid is None:
        grid = RasterModelGrid(
            shape, xy_spacing=xy_spacing, xy_of_lower_left=xy_of_lower_left
        )
    if name:
        grid.add_field(name, data, at="node")

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
    >>> import os
    >>> from landlab import RasterModelGrid
    >>> from landlab.io.esri_ascii import write_esri_ascii

    >>> grid = RasterModelGrid((4, 5), xy_spacing=(2.0, 2.0))
    >>> grid.at_node["air__temperature"] = np.arange(20.0)
    >>> files = write_esri_ascii("test.asc", grid)  # doctest: +SKIP
    >>> [os.path.basename(name) for name in sorted(files)]  # doctest: +SKIP
    ['test.asc']

    >>> _ = grid.add_field("land_surface__elevation", np.arange(20.0), at="node")
    >>> grid.at_node["land_surface__elevation"] = np.arange(20.0)
    >>> files = write_esri_ascii("test.asc", grid)  # doctest: +SKIP
    >>> [os.path.basename(name) for name in sorted(files)]  # doctest: +SKIP
    ['test_air__temperature.asc', 'test_land_surface__elevation.asc']
    """
    if os.path.exists(path) and not clobber:
        raise ValueError("file exists")

    if isinstance(names, (str, pathlib.Path)):
        names = [names]

    names = names or fields.at_node.keys()

    if len(names) == 1:
        paths = [path]
    elif len(names) > 1:
        (base, ext) = os.path.splitext(path)
        paths = [base + "_" + name + ext for name in names]
    else:
        raise ValueError("no node fields to write")

    bad_names = set(names) - set(fields.at_node.keys())
    if len(bad_names) > 0:
        raise ValueError("unknown field name(s): %s" % ",".join(bad_names))

    header = {
        "ncols": fields.number_of_node_columns,
        "nrows": fields.number_of_node_rows,
        "xllcorner": fields.node_x[0],
        "yllcorner": fields.node_y[0],
        "cellsize": fields.dx,
    }

    for path, name in zip(paths, names):
        header_lines = [f"{key} {str(val)}" for key, val in list(header.items())]
        data = fields.at_node[name].reshape(header["nrows"], header["ncols"])
        np.savetxt(
            path, np.flipud(data), header=os.linesep.join(header_lines), comments=""
        )

    return paths

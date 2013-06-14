#! /usr/bin/env python
"""
Read data from an ESRI ASCII file into a RasterModelGrid
"""

import types

import numpy as np

from landlab import RasterModelGrid


def _read_asc_header(asc_file):
    """
    Read header information from an ESRI ASCII raster file.

    The header contains the following variables,
        :ncols: Number of columns
        :nrows: Number of rows
        :xllcorner: X (column) coordinate of lower-left coordinate of grid
        :yllcorner: Y (row) coordinate of lower-left coordinate of grid
        :cellsize: Grid spacing between rows and columns
        :nodata_value: No-data value

    .. todo::
        - nodata_value is optional
        - xllcorner or xllcenter
        - yllcorner or yllcenter
    """
    from itertools import islice

    required_keys = set(['ncols', 'nrows', 'xllcorner', 'yllcorner',
                         'cellsize', 'NODATA_value'])

    header = dict()
    for line in list(islice(asc_file, len(required_keys))):
        items = line.split()
        if items[0] in required_keys:
            header[items[0]] = items[1]
        else:
            raise KeyError(items[0])

    return header


def _read_asc_data(asc_file):
    """
    Read gridded data from an ESRI ASCII data file.

    :asc_file: File-like object of the data file pointing to the start of the
               data.

    .. note::
        Row 1 of the data is at the top of the raster, row 2 is just under
        row 1, and so on.
    """
    return np.genfromtxt(asc_file)


def read_esri_ascii(asc_file):
    """
    Read data from an ESRI_ ASCII file into a RasterModelGrid.

    :asc_file: Name of the ESRI ASCII file of a file-like object

    .. _ESRI: http://resources.esri.com/help/9.3/arcgisengine/java/GP_ToolRef/spatial_analyst_tools/esri_ascii_raster_format.htm
    """
    if isinstance(asc_file, types.StringTypes):
        file_name = asc_file
        with open(file_name, 'r') as asc_file:
            header = _read_asc_header(asc_file)
            data = _read_asc_data(asc_file)
    else:
        file_name = asc_file.name
        header = _read_asc_header(asc_file)
        data = _read_asc_data(asc_file)

    shape = (int(header['nrows']), int(header['ncols']))
    spacing = (float(header['cellsize']), float(header['cellsize']))
    origin = (float(header['xllcorner']), float(header['yllcorner']))

    grid = RasterModelGrid(num_rows=shape[0], num_cols=shape[1], dx=spacing[0])

    return (grid, data)

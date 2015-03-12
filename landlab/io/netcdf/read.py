#! /usr/bin/env python
"""
Read data from a NetCDF file into a RasterModelGrid.
"""

try:
    import netCDF4 as nc4
except ImportError:
    import warnings
    warnings.warn('Unable to import netCDF4.', ImportWarning)

from scipy.io import netcdf as nc

import os
import types
import re

import numpy as np

from landlab.io.netcdf.errors import NotRasterGridError
from landlab.io.netcdf._constants import (_AXIS_DIMENSION_NAMES,
                                          _AXIS_COORDINATE_NAMES,
                                          _COORDINATE_NAMES)


def _length_of_axis_dimension(root, axis_name):
    try:
        return len(root.dimensions[axis_name])
    except TypeError:
        return root.dimensions[axis_name]


def _read_netcdf_grid_shape(root):
    shape = []
    for axis_name in _AXIS_DIMENSION_NAMES:
        try:
            shape.append(_length_of_axis_dimension(root, axis_name))
        except KeyError:
            pass
    return shape


def _read_netcdf_coordinate_values(root):
    values = []
    for coordinate_name in _AXIS_COORDINATE_NAMES:
        try:
            values.append(root.variables[coordinate_name][:].copy())
        except KeyError:
            pass
    return values


def _read_netcdf_coordinate_units(root):
    units = []
    for coordinate_name in _AXIS_COORDINATE_NAMES:
        try:
            units.append(root.variables[coordinate_name].units)
        except KeyError:
            pass
    return units


def _read_netcdf_structured_grid(root):
    shape = _read_netcdf_grid_shape(root)
    coordinates = _read_netcdf_coordinate_values(root)
    units = _read_netcdf_coordinate_units(root)

    for coordinate in coordinates:
        coordinate.shape = shape

    return coordinates


def _read_netcdf_structured_data(root):
    fields = dict()
    for (name, var) in root.variables.items():
        if not name in _COORDINATE_NAMES:
            fields[name] = var[:].copy()
            fields[name].shape = (fields[name].size, )
    return fields


def _get_raster_spacing(coords):
    spacing = np.empty(len(coords), dtype=np.float64)

    for (axis, coord) in enumerate(coords):
        coord_spacing = np.diff(coord, axis=axis)
        try:
            assert(np.all(coord_spacing == coord_spacing.flat[0]))
        except AssertionError:
            raise NotRasterGridError()
        spacing[axis] = coord_spacing.flat[0]

    try:
        assert(np.all(spacing == spacing[0]))
    except AssertionError:
        raise NotRasterGridError()
    else:
        return spacing[0]


def read_netcdf(nc_file, reshape=False, just_grid=False):
    """
    Reads the NetCDF file *nc_file*, and writes it to the fields of a new
    RasterModelGrid, which it then returns.
    Check the names of the fields in the returned grid with
    grid.at_nodes.keys().
    """
    from landlab import RasterModelGrid

    try:
        root = nc.netcdf_file(nc_file, 'r', version=2)
    except TypeError:
        root = nc4.Dataset(nc_file, 'r', format='NETCDF4')

    node_coords = _read_netcdf_structured_grid(root)

    assert(len(node_coords) == 2)

    spacing = _get_raster_spacing(node_coords)

    shape = node_coords[0].shape

    grid = RasterModelGrid(num_rows=shape[0], num_cols=shape[1],
                           dx=spacing)

    if not just_grid:
        fields = _read_netcdf_structured_data(root)
        for (name, values) in fields.items():
            grid.add_field('node', name, values)

    root.close()

    return grid

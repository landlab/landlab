#! /usr/bin/env python
"""
Read data from a GEBCO NetCDF file into a RasterModelGrid.
"""

try:
    import netCDF4 as nc4
except ImportError:
    import warnings
    warnings.warn('Unable to import netCDF4.', ImportWarning)

from scipy.io import netcdf as nc

import numpy as np

from landlab import RasterModelGrid
from landlab.io.netcdf.errors import NotRasterGridError


_COORDINATE_NAMES = ['x_range', 'y_range', 'z_range', 'spacing', 'dimension']
_AXIS_NAMES = ['x', 'y']


def _read_netcdf_grid_shape(root):
    return root.variables['dimension'][:]


def _read_netcdf_grid_spacing(root):
    return root.variables['spacing'][:]


def _read_netcdf_coordinate_values(root):
    values = []

    spacing = root.variables['spacing'][:]
    for (axis, coordinate_name) in enumerate(_AXIS_NAMES):
        coord_range = root.variables[coordinate_name + '_range'][:]
        values.append(np.arange(coord_range[0], coord_range[1] + spacing[axis],
                                spacing[axis]))

    return values


def _read_netcdf_coordinate_units(root):
    return ['degrees_east', 'degrees_north']


def _read_netcdf_structured_data(root):
    fields = dict()
    for (name, var) in root.variables.items():
        if name not in _COORDINATE_NAMES:
            fields[name] = var[:] * var.scale_factor + var.add_offset
            fields[name].shape = (fields[name].size, )
    return fields


def read_netcdf(nc_file, reshape=False, just_grid=False):
    """
    Reads the NetCDF file *nc_file*, and writes it to the fields of a new
    RasterModelGrid, which it then returns.
    Check the names of the fields in the returned grid with
    grid.at_nodes.keys().
    """
    try:
        root = nc.netcdf_file(nc_file, 'r', version=2)
    except TypeError:
        root = nc4.Dataset(nc_file, 'r', format='NETCDF4')

    shape = _read_netcdf_grid_shape(root)
    spacing = _read_netcdf_grid_spacing(root)

    assert(len(shape) == 2)
    assert(len(spacing) == 2)
    if spacing[0] != spacing[1]:
        raise NotRasterGridError()

    grid = RasterModelGrid(num_rows=shape[0], num_cols=shape[1], dx=spacing[0])

    if not just_grid:
        fields = _read_netcdf_structured_data(root)
        for (name, values) in fields.items():
            grid.add_field('node', name, values)

    root.close()

    return grid

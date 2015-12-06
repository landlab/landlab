#! /usr/bin/env python
"""Read data from a GEBCO NetCDF file into a RasterModelGrid."""


try:
    import netCDF4 as nc4
except ImportError:
    import warnings
    warnings.warn('Unable to import netCDF4.', ImportWarning)

from scipy.io import netcdf as nc

from landlab import RasterModelGrid
from landlab.io.netcdf.errors import NotRasterGridError


_COORDINATE_NAMES = ['x_range', 'y_range', 'z_range', 'spacing', 'dimension']
_AXIS_NAMES = ['x', 'y']


def _read_netcdf_grid_shape(root):
    """Read the grid shape from a GEBCO NetCDF file.

    Parameters
    ----------
    root : netcdf_file
        A NetCDF file.

    Returns
    -------
    tuple of int
        The shape of the grid as number of rows, then columns.
    """
    return root.variables['dimension'][:]


def _read_netcdf_grid_spacing(root):
    """Read the grid spacing from a GEBCO NetCDF file.

    Parameters
    ----------
    root : netcdf_file
        A NetCDF file.

    Returns
    -------
    tuple of float
        The spacing of the grid between rows, then columns.
    """
    return root.variables['spacing'][:]


def _read_netcdf_structured_data(root):
    """Read the grid data from a GEBCO NetCDF file.

    Parameters
    ----------
    root : netcdf_file
        A NetCDF file.

    Returns
    -------
    dict
        The data fields as numpy arrays. Keys are the variable names, and
        values are the data.
    """
    fields = dict()
    for (name, var) in root.variables.items():
        if name not in _COORDINATE_NAMES:
            fields[name] = var[:] * var.scale_factor + var.add_offset
            fields[name].shape = (fields[name].size, )
    return fields


def read_netcdf(nc_file, just_grid=False):
    """Read a GEBCO-formatted NetCDF file.

    Reads the NetCDF file *nc_file*, and writes it to the fields of a new
    RasterModelGrid, which it then returns.  Check the names of the fields
    in the returned grid with grid.at_nodes.keys().

    Parameters
    ----------
    nc_file : str
        Name of the NetCDF file.
    just_grid : bool, optional
        If ``True``, just read the grid information and forget the data.
        Otherwise add the data as fields.

    Returns
    -------
    RasterModelGrid
        A newly-created :any:`RasterModelGrid`.
    """
    try:
        root = nc.netcdf_file(nc_file, 'r', version=2)
    except TypeError:
        root = nc4.Dataset(nc_file, 'r', format='NETCDF4')

    shape = _read_netcdf_grid_shape(root)
    spacing = _read_netcdf_grid_spacing(root)

    assert len(shape) == 2
    assert len(spacing) == 2
    if spacing[0] != spacing[1]:
        raise NotRasterGridError()

    grid = RasterModelGrid(shape, spacing=spacing)

    if not just_grid:
        fields = _read_netcdf_structured_data(root)
        for (name, values) in fields.items():
            grid.add_field('node', name, values)

    root.close()

    return grid

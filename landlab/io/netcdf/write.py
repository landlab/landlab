#! /usr/bin/env python
import os
import warnings
import types


try:
    import netCDF4 as nc4
except ImportError:
    warnings.warn('Unable to import netCDF4.', ImportWarning)

from scipy.io import netcdf as nc


from landlab.io.netcdf._constants import (_AXIS_DIMENSION_NAMES,
                                          _AXIS_COORDINATE_NAMES,
                                          _NP_TO_NC_TYPE)

def _set_netcdf_attributes(root, attrs):
    """
    Set attributes of the netCDF Database object, *root*. Attributes are
    given as key/value pairs from *attrs*.
    """
    for (key, val) in attrs.items():
        setattr(root, key, val)


def _get_dimension_names(shape):
    names = _AXIS_DIMENSION_NAMES[- 1: - (len(shape) + 1): - 1]
    return names[::-1]


def _get_dimension_sizes(shape):
    names = _AXIS_DIMENSION_NAMES[- 1: - (len(shape) + 1): - 1]

    sizes = dict()
    for (axis, name) in enumerate(names):
        sizes[name] = shape[- (axis + 1)]

    return sizes


def _get_axes_names(shape):
    names = _AXIS_COORDINATE_NAMES[- 1: - (len(shape) + 1): - 1]
    return names[::-1]


def _set_netcdf_structured_dimensions(root, shape):
    """
    Add dimensions to *root* for a structured grid of size *shape*. The
    dimension names will be 'ni', 'nj', and 'nk'. 'ni' is the length of the
    fast dimension, followed by 'nj', and then 'nk'.

    For example, a grid with shape (3, 4, 5) will have dimensions ni=5,
    nj=4, and nk=3. Lower dimension grids simply drop the slowest dimension.
    Thus, a grid with shape (3, 4) has dimensions ni=4, and nj=3.
    """
    assert(len(shape) >= 1 and len(shape) <= 3)

    dimensions = _get_dimension_sizes(shape)

    dims = root.dimensions
    for (name, dim_size) in dimensions.items():
        if name not in dims:
            root.createDimension(name, dim_size)

    if not 'nt' in dims:
        nt = root.createDimension('nt', None)


def _set_netcdf_variables(root, fields, **kwds):
    """Set the field variables.
    
    First set the variables that define the grid and then the variables at
    the grid nodes and cells.
    """
    names = kwds.pop('names', None)

    _add_spatial_variables(root, fields, **kwds)
    #_add_variables_at_points(root, fields['nodes'])
    _add_variables_at_points(root, fields, names=names)


def _add_spatial_variables(root, grid, **kwds):
    """Add the variables to *root* that define the structured grid, *grid*.
    """
    long_name = kwds.get('long_name', {})

    vars = root.variables
    dims = root.dimensions

    spatial_variable_names = _get_axes_names(grid.shape)
    spatial_variable_shape = _get_dimension_names(grid.shape)

    for (axis, name) in enumerate(spatial_variable_names):
        try:
            var = vars[name]
        except KeyError:
            var = root.createVariable(name, 'f8', spatial_variable_shape)

        coords = grid.node_axis_coordinates(axis=axis).view()
        coords.shape = var.shape
        var[:] = coords

        var.units = grid.axis_units[axis]
        try:
            var.long_name = long_name[name]
        except KeyError:
            var.long_name = grid.axis_name[axis]


def _add_variables_at_points(root, fields, names=None):
    if isinstance(names, types.StringTypes):
        names = [names]
    names = names or fields['node'].keys()

    vars = root.variables

    spatial_variable_shape = _get_dimension_names(fields.shape)

    try:
        n_times = len(vars['t']) - 1
    except KeyError:
        n_times = 0

    node_fields = fields['node']
    for var_name in names:
        try:
            var = vars[var_name]
        except KeyError:
            var = root.createVariable(
                var_name, _NP_TO_NC_TYPE[str(node_fields[var_name].dtype)],
                ['nt'] + spatial_variable_shape)

        if node_fields[var_name].size > 1:
            data = node_fields[var_name].view()
            data.shape = var.shape[1:]
            try:
                var[n_times, :] = data
            except ValueError:
                raise
        else:
            var[n_times] = node_fields[var_name].flat[0]

        var.units = node_fields.units[var_name] or '?'
        var.long_name = var_name


def _add_time_variable(root, time, **kwds):
    units = kwds.get('units', 'days')
    reference = kwds.get('reference', '00:00:00 UTC')

    vars = root.variables

    try:
        t = vars['t']
    except KeyError:
        t = root.createVariable('t', 'f8', ('nt', ))
        t.units = ' '.join([units, 'since', reference])
        t.long_name = 'time'

    n_times = len(t)
    if time is not None:
        t[n_times] = time
    else:
        t[n_times] = n_times


_VALID_NETCDF_FORMATS = set([
    'NETCDF3_CLASSIC',
    'NETCDF3_64BIT',
    'NETCDF4_CLASSIC',
    'NETCDF4',
])

def write_netcdf(path, fields, attrs=None, append=False,
                 format='NETCDF3_64BIT', names=None):
    """Write landlab fields to netcdf.

    Write the data and grid information for *fields* to *path* as NetCDF.
    If the *append* keyword argument in True, append the data to an existing
    file, if it exists. Otherwise, clobber an existing files.

    Parameters
    ----------
    path : str
        Path to output file.
    fields : field-like
        Landlab field object that holds a grid and associated values.
    append : boolean, optional
        Append data to an existing file, otherwise clobber the file.
    format : {'NETCDF3_CLASSIC', 'NETCDF3_64BIT', 'NETCDF4_CLASSIC', 'NETCDF4'}
        Format of output netcdf file.
    attrs : dict
        Attributes to add to netcdf file.
    names : iterable of str, optional
        Names of the fields to include in the netcdf file. If not provided,
        write all fields.

    Examples
    --------
    >>> from landlab import RasterModelGrid
    >>> from landlab.io.netcdf import write_netcdf

    Create a uniform rectilinear grid with four rows and 3 columns, and add
    some data fields to it.

    >>> rmg = RasterModelGrid(4, 3)
    >>> _ = rmg.add_field('node', 'topographic_elevation', np.arange(12.))
    >>> _ = rmg.add_field('node', 'uplift_rate', 2. * np.arange(12.))

    Create a temporary directory to write the netcdf file into.

    >>> import tempfile, os
    >>> temp_dir = tempfile.mkdtemp()
    >>> os.chdir(temp_dir)

    Write the grid to a netcdf3 file but only include the *uplift_rate*
    data in the file.

    >>> write_netcdf('test.nc', rmg, format='NETCDF3_64BIT',
    ...     names='uplift_rate')

    Read the file back in and check its contents.

    >>> from scipy.io import netcdf
    >>> fp = netcdf.netcdf_file('test.nc', 'r')
    >>> 'uplift_rate' in fp.variables
    True
    >>> 'topographic_elevation' in fp.variables
    False
    >>> fp.variables['uplift_rate'][:].flatten()
    array([  0.,   2.,   4.,   6.,   8.,  10.,  12.,  14.,  16.,  18.,  20.,
           22.])
    """
    if format not in _VALID_NETCDF_FORMATS:
        raise ValueError('format not understood')

    attrs = attrs or {}

    if os.path.isfile(path) and append:
        mode = 'a'
    else:
        mode = 'w'

    if format == 'NETCDF3_CLASSIC':
        root = nc.netcdf_file(path, mode, version=1)
    elif format == 'NETCDF3_64BIT':
        root = nc.netcdf_file(path, mode, version=2)
    else:
        root = nc4.Dataset(path, mode, format=format)

    _set_netcdf_attributes(root, attrs)
    _set_netcdf_structured_dimensions(root, fields.shape)
    _set_netcdf_variables(root, fields, names=names)

    root.close()

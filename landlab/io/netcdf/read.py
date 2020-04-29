#! /usr/bin/env python
"""Read data from a NetCDF file into a RasterModelGrid.

Read netcdf
+++++++++++

.. autosummary::

    ~landlab.io.netcdf.read.read_netcdf
"""
import numpy as np
import xarray as xr

from landlab.io import (
    MismatchGridDataSizeError,
    MismatchGridXYLowerLeft,
    MismatchGridXYSpacing,
)
from landlab.io.netcdf._constants import (
    _AXIS_COORDINATE_NAMES,
    _AXIS_DIMENSION_NAMES,
    _COORDINATE_NAMES,
)
from landlab.io.netcdf.errors import NotRasterGridError
from landlab.utils import add_halo


def _length_of_axis_dimension(root, axis_name):
    """Get the size of an axis by axis name.

    Parameters
    ----------
    root : netcdf_file
        A NetCDF object.
    axis_name : str
        Name of the axis in the NetCDF file.

    Returns
    -------
    int
        Size of the dimension.
    """
    try:
        return len(root.dimensions[axis_name])
    except TypeError:
        return root.dimensions[axis_name]


def _read_netcdf_grid_shape(root):
    """Get the shape of a raster grid.

    Parameters
    ----------
    root : netcdf_file
        A NetCDF object.

    Returns
    -------
    (rows, columns)
        Shape of the grid.
    """
    shape = []
    for axis_name in _AXIS_DIMENSION_NAMES:
        try:
            shape.append(_length_of_axis_dimension(root, axis_name))
        except KeyError:
            pass
    return tuple(shape)


def _read_netcdf_coordinate_values(root):
    """Get arrays of coordinates for grid points.

    Parameters
    ----------
    root : netcdf_file
        A NetCDF file.

    Returns
    -------
    tuple of ndarray
        Node coordinates for each dimension.
    """
    values = []
    for coordinate_name in _AXIS_COORDINATE_NAMES:
        try:
            values.append(root.variables[coordinate_name][:].copy())
        except KeyError:
            pass
    return tuple(values)


def _read_netcdf_coordinate_units(root):
    """Get units for coodinate values.

    Parameters
    ----------
    root : netcdf_file
        A NetCDF file.

    Returns
    -------
    tuple of str
        Units for each coordinate.
    """
    units = []
    for coordinate_name in _AXIS_COORDINATE_NAMES:
        try:
            units.append(root.variables[coordinate_name].units)
        except KeyError:
            pass
    return tuple(units)


def _read_netcdf_structured_grid(root):
    """Get node coordinates for a structured grid.

    Parameters
    ----------
    root : netcdf_file
        A NetCDF file.

    Returns
    -------
    tuple of ndarray
        Node coordinates for each dimension reshaped to match the shape
        of the grid.
    """
    shape = _read_netcdf_grid_shape(root)
    coordinates = _read_netcdf_coordinate_values(root)

    for coordinate in coordinates:
        coordinate.shape = shape

    return coordinates


def _read_netcdf_raster_structured_grid(root):
    """Get node coordinates for a structured grid written as a raster.

    Parameters
    ----------
    root : netcdf_file
        A NetCDF file.

    Returns
    -------
    tuple of ndarray
        Node coordinates for each dimension reshaped to match the shape
        of the grid.
    """
    shape = _read_netcdf_grid_shape(root)
    coordinates = _read_netcdf_coordinate_values(root)

    if len(coordinates) != 2:
        assert ValueError("Rasters must have only two spatial coordinate dimensions")
    else:
        coordinates = np.meshgrid(coordinates[0], coordinates[1], indexing="ij")

    for coordinate in coordinates:
        coordinate.shape = shape

    return coordinates


def _read_netcdf_structured_data(root):
    """Get data values for a structured grid.

    Parameters
    ----------
    root : netcdf_file
        A NetCDF file.

    Returns
    -------
    dict
        Data values, reshaped to match that of the grid. Keys are the
        variable names as read from the NetCDF file.
    """
    fields = dict()
    grid_mapping_exists = False
    grid_mapping_dict = None
    for (name, var) in root.variables.items():
        # identify if a grid mapping variable exist and do not pass it as a field
        if name not in _COORDINATE_NAMES:
            if hasattr(var, "grid_mapping"):
                grid_mapping = getattr(var, "grid_mapping")
                if type(grid_mapping) is bytes:
                    grid_mapping = grid_mapping.decode("utf-8")
                grid_mapping_exists = True

    dont_use = list(_COORDINATE_NAMES)
    if grid_mapping_exists:
        dont_use.append(grid_mapping)
    for (name, var) in root.variables.items():
        if name not in dont_use:
            fields[name] = var.values.reshape((-1,))

    if grid_mapping_exists:
        grid_mapping_variable = root.variables[grid_mapping]
        grid_mapping_dict = {"name": grid_mapping}
        try:
            for att in grid_mapping_variable.ncattrs():
                grid_mapping_dict[att] = getattr(grid_mapping_variable, att)
        except AttributeError:  # if scipy is doing the reading
            for att in var._attributes:
                grid_mapping_dict[att] = getattr(grid_mapping_variable, att)
    return fields, grid_mapping_dict


def _get_raster_spacing(coords):
    """Get the row and column spacing of a raster grid.

    Parameters
    ----------
    root : netcdf_file
        A NetCDF object.

    Returns
    -------
    (dy, dx)
        Spacing of grid rows and columns.
    """
    spacing = np.empty(len(coords), dtype=np.float64)

    for (axis, coord) in enumerate(coords):
        coord_spacing = np.diff(coord, axis=axis)
        if not np.all(coord_spacing == coord_spacing.flat[0]):
            raise NotRasterGridError()
        spacing[axis] = coord_spacing.flat[0]

    if not np.all(spacing == spacing[0]):
        raise NotRasterGridError()

    return spacing[0]


def read_netcdf(
    nc_file, grid=None, name=None, just_grid=False, halo=0, nodata_value=-9999.0,
):
    """Create a :class:`~.RasterModelGrid` from a netcdf file.

    Create a new :class:`~.RasterModelGrid` from the netcdf file, *nc_file*.
    If the netcdf file also contains data, add that data to the grid's fields.
    To create a new grid without any associated data from the netcdf file,
    set the *just_grid* keyword to ``True``.

    A halo can be added with the keyword *halo*.

    If you want the fields to be added to an existing grid, it can be passed
    to the keyword argument *grid*.

    Parameters
    ----------
    nc_file : str
        Name of a netcdf file.
    grid : *grid* , optional
        Adds data to an existing *grid* instead of creating a new one.
    name : str, optional
        Add only fields with NetCDF variable name to the grid. Default is to
        add all NetCDF varibles to the grid.
    just_grid : boolean, optional
        Create a new grid but don't add value data.
    halo : integer, optional
        Adds outer border of depth halo to the *grid*.
    nodata_value : float, optional
        Value that indicates an invalid value. Default is -9999.

    Returns
    -------
    :class:`~.RasterModelGrid`
        A newly-created :class:`~.RasterModelGrid`.

    Examples
    --------
    Import :func:`read_netcdf` and the path to an example netcdf file included
    with landlab.

    >>> from landlab.io.netcdf import read_netcdf

    Create a new grid from the netcdf file. The example grid is a uniform
    rectilinear grid with 4 rows and 3 columns of nodes with unit spacing.
    The data file also contains data defined at the nodes for the grid for
    a variable called, *surface__elevation*.

    >>> grid = read_netcdf("test-netcdf4.nc") # doctest: +SKIP
    >>> grid.shape == (4, 3) # doctest: +SKIP
    True
    >>> grid.dy, grid.dx # doctest: +SKIP
    (1.0, 1.0)
    >>> list(grid.at_node.keys()) # doctest: +SKIP
    ['surface__elevation']
    >>> grid.at_node['surface__elevation'] # doctest: +SKIP
    array([  0.,   1.,   2.,   3.,   4.,   5.,   6.,   7.,   8.,   9.,  10.,
            11.])

    :func:`read_netcdf` will try to determine the format of the netcdf file.
    For example, the same call will also work for *netcdf3*-formatted files.

    >>> grid = read_netcdf("test-netcdf3-64bit.nc") # doctest: +SKIP
    >>> grid.shape == (4, 3) # doctest: +SKIP
    True
    >>> grid.dy, grid.dx # doctest: +SKIP
    (1.0, 1.0)

    A more complicated example might add data with a halo to an existing grid.
    Note that the lower left corner must be specified correctly for the data
    and the grid to align correctly.

    >>> from landlab import RasterModelGrid
    >>> grid = RasterModelGrid((6, 5), xy_of_lower_left=(-1., -1.)) # doctest: +SKIP
    >>> grid = read_netcdf(
    ...     "test-netcdf4.nc",
    ...     grid=grid,
    ...     halo=1,
    ...     nodata_value=-1,
    ... ) # doctest: +SKIP
    >>> grid.at_node['surface__elevation'].reshape(grid.shape) # doctest: +SKIP
    array([[ -1.,  -1.,  -1.,  -1.,  -1.],
           [ -1.,   0.,   1.,   2.,  -1.],
           [ -1.,   3.,   4.,   5.,  -1.],
           [ -1.,   6.,   7.,   8.,  -1.],
           [ -1.,   9.,  10.,  11.,  -1.],
           [ -1.,  -1.,  -1.,  -1.,  -1.]])
    """
    from landlab import RasterModelGrid

    dataset = xr.open_dataset(nc_file)

    if isinstance(name, str):
        names = {name}
    elif name is None:
        names = set(dataset.variables)
    else:
        names = set(name)

    # test if the input is a raster (x and y) are only 1-D instead of 2D.
    if len(dataset["x"].shape) == 1:
        y, x = np.meshgrid(dataset["y"], dataset["x"], indexing="ij")
    else:
        x = dataset["x"]
        y = dataset["y"]

    dx = np.diff(x, axis=1)
    dy = np.diff(y, axis=0)

    if np.all(dx == dx[0, 0]) and np.all(dy == dy[0, 0]):
        xy_spacing = (dx[0, 0], dy[0, 0])
    else:
        raise NotRasterGridError()

    shape = x.shape
    xy_of_lower_left = (
        x[0, 0] - halo * xy_spacing[0],
        y[0, 0] - halo * xy_spacing[1],
    )

    if grid is None:
        grid = RasterModelGrid(
            shape, xy_spacing=xy_spacing, xy_of_lower_left=xy_of_lower_left
        )
    else:
        if grid.shape != (shape[0] + 2 * halo, shape[1] + 2 * halo):
            raise MismatchGridDataSizeError(
                shape[0] + 2 * halo * shape[1] + 2 * halo,
                grid.number_of_node_rows * grid.number_of_node_columns,
            )
        if (grid.dx, grid.dy) != xy_spacing:
            raise MismatchGridXYSpacing((grid.dx, grid.dy), xy_spacing)

        if grid.xy_of_lower_left != xy_of_lower_left:
            raise MismatchGridXYLowerLeft(grid.xy_of_lower_left, xy_of_lower_left)

    if not just_grid:
        fields, grid_mapping_dict = _read_netcdf_structured_data(dataset)
        for (field_name, values) in fields.items():

            # add halo if necessary
            if halo > 0:
                values = add_halo(
                    values.reshape(shape), halo=halo, halo_value=nodata_value
                ).reshape((-1,))

            # add only the requested fields.
            if (name is None) or (field_name == name):
                add_field = True
            else:
                add_field = False

            if add_field:
                grid.add_field(field_name, values, at="node", clobber=True)

        if (name is not None) and (name not in grid.at_node):
            raise ValueError(
                "Specified field {name} was not in provided NetCDF.".format(name=name)
            )

    ignore = {"x", "y"}
    for name in names - ignore:
        values = dataset.variables[name].values
        if halo > 0:
            values = add_halo(
                values.reshape(shape), halo=halo, halo_value=nodata_value
            ).reshape((-1,))
        grid.add_field(name, values, at="node", clobber=True)

    return grid

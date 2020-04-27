import fnmatch

import xarray as xr

from ...grid.hex import HexModelGrid
from ...grid.radial import RadialModelGrid
from ...grid.raster import RasterModelGrid


def from_netcdf(filename_or_obj, include="*", exclude=None):
    """Create a :class:`~.ModelGrid` from a netcdf file.

    Create a new :class:`~.ModelGrid` from the netcdf file, *nc_file*.
    If the netcdf file also contains data, add that data to the grid's fields.
    To create a new grid without any associated data from the netcdf file,
    use *include=None*.

    Parameters
    ----------
    filename_or_obj : str, Path, or file
        Strings and Path objects are interpreted as a path to a netCDF file
        or an OpenDAP URL and opened with python-netCDF4, unless the
        filename ends with .gz, in which case the file is gunzipped and
        opened with scipy.io.netcdf (only netCDF3 supported). Byte-strings
        or file-like objects are opened by scipy.io.netcdf (netCDF3)
        or h5py (netCDF4/HDF).
    include : str or iterable of str, optional
        A list of unix-style glob patterns of field names to include. Fully
        qualified field names that match any of these patterns will be
        written to the output file. A fully qualified field name is one that
        that has a prefix that indicates what grid element is defined on
        (e.g. "at_node:topographic__elevation"). The default is to include
        all fields.
    exclude : str or iterable of str, optional
        Like the *include* keyword but, instead, fields matching these
        patterns will be excluded from the output file.

    Returns
    -------
    :class:`~.ModelGrid`
        A newly-created :class:`~.ModelGrid`.
    """
    include = include or []
    if isinstance(include, str):
        include = [include]
    exclude = exclude or []
    if isinstance(exclude, str):
        exclude = [exclude]

    _grid_from_str = {
        "RasterModelGrid": RasterModelGrid,
        "HexModelGrid": HexModelGrid,
        "RadialModelGrid": RadialModelGrid,
        "uniform_rectilinear": RasterModelGrid,
        "triangular": HexModelGrid,
        "radial": RadialModelGrid,
    }

    with xr.open_dataset(filename_or_obj) as dataset:
        grid_type = dataset.attrs["grid_type"]

        try:
            from_dataset = _grid_from_str[grid_type].from_dataset
        except KeyError:
            raise RuntimeError("grid type not recognized ({0})".format(grid_type))
        else:
            grid = from_dataset(dataset)

        qualified_names = [name for name in dataset.variables if name.startswith("at_")]
        names = set()
        for pattern in include:
            names.update(fnmatch.filter(qualified_names, pattern))
        for pattern in exclude:
            names.difference_update(fnmatch.filter(qualified_names, pattern))

        for name in names:
            at_name, field_name = name.split(":")
            getattr(grid, at_name)[field_name] = dataset[name]

        grid.status_at_node = dataset["status_at_node"]

    return grid

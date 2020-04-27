import pathlib

import numpy as np
import xarray as xr


def to_netcdf(
    grid, path, include="*", exclude=None, time=None, format="NETCDF4", mode="w"
):
    """Write landlab a grid to a netcdf file.

    Write the data and grid information for *grid* to *path* as NetCDF.
    If the *append* keyword argument in True, append the data to an existing
    file, if it exists. Otherwise, clobber an existing files.

    Parameters
    ----------
    grid : ModelGrid
        Landlab grid object that holds a grid and field values.
    path : str
        Path to which to save this grid.
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
    format : {'NETCDF3_CLASSIC', 'NETCDF3_64BIT', 'NETCDF4_CLASSIC', 'NETCDF4'}
        Format of output netcdf file.
    attrs : dict
        Attributes to add to netcdf file.
    mode : {"w", "a"}, optional
        Write ("w") or append ("a") mode. If mode="w", any existing file at
        this location will be overwritten. If mode="a", existing variables
        will be overwritten.


    Parameters
    ----------
    Examples
    --------
    >>> import numpy as np
    >>> from landlab import RasterModelGrid
    >>> from landlab.io.netcdf import to_netcdf

    Create a uniform rectilinear grid with four rows and 3 columns, and add
    some data fields to it.

    >>> rmg = RasterModelGrid((4, 3))
    >>> rmg.at_node["topographic__elevation"] = np.arange(12.0)
    >>> rmg.at_node["uplift_rate"] = 2.0 * np.arange(12.0)

    Create a temporary directory to write the netcdf file into.

    >>> import tempfile, os
    >>> temp_dir = tempfile.mkdtemp()
    >>> os.chdir(temp_dir)

    Write the grid to a netcdf3 file but only include the *uplift_rate*
    data in the file.

    >>> to_netcdf(
    ...     rmg, "test.nc", format="NETCDF3_64BIT", include="at_node:uplift_rate"
    ... )

    Read the file back in and check its contents.

    >>> from scipy.io import netcdf
    >>> fp = netcdf.netcdf_file('test.nc', 'r')
    >>> 'at_node:uplift_rate' in fp.variables
    True
    >>> 'at_node:topographic__elevation' in fp.variables
    False
    >>> fp.variables['at_node:uplift_rate'][:].flatten()
    array([  0.,   2.,   4.,   6.,   8.,  10.,  12.,  14.,  16.,  18.,  20.,
            22.])

    >>> rmg.at_cell["air__temperature"] = np.arange(2.0)
    >>> to_netcdf(
    ...     rmg,
    ...     "test-cell.nc",
    ...     format="NETCDF3_64BIT",
    ...     include="at_cell:*",
    ...     # names="air__temperature", at="cell",
    ... )
    """
    path = pathlib.Path(path)

    ds = grid.as_dataset(include=include, exclude=exclude)
    new_fields = set([name for name in ds.variables if name.startswith("at_")])
    if time is not None:
        for name in new_fields:
            dim, _ = name[len("at_") :].split(":")
            ds[name] = (("time",) + ds[name].dims, ds[name].values[None])
        ds["time"] = (("time",), [time])

    if format != "NETCDF4":
        ds["status_at_node"] = (
            ("node",),
            ds["status_at_node"].values.astype(dtype=int),
        )

    if mode == "a" and path.is_file():
        with xr.open_dataset(path) as dataset:
            try:
                times = dataset["time"].values
            except KeyError:
                times = [0.0]
            if time is None:
                time = times[-1] + 1.0

            data = {
                "time": (("time",), np.concatenate([times, [time]])),
            }

            existing_fields = set(
                [name for name in dataset.variables if name.startswith("at_")]
            )
            for name in existing_fields:
                dims, array = dataset[name].dims, dataset[name].values
                if "time" not in dims:
                    dims, array = ("time",) + dims, array[None]
                data[name] = (dims, np.concatenate([array, ds[name].values]))

            for name in new_fields - existing_fields:
                dims, array = ds[name].dims, ds[name].values
                if "time" not in dims:
                    dims, array = ("time",) + dims, array[None]
                padded_array = np.empty_like(
                    array, shape=(array.shape[0] + 1,) + array.shape[1:],
                )
                padded_array[-1] = array
                data[name] = (dims, padded_array)

            ds.update(data)

    ds.to_netcdf(path, format=format, mode="w", unlimited_dims=("time",))

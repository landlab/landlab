import os

from .errors import NotRasterGridError
from .read import read_netcdf
from .write import write_netcdf, write_raster_netcdf

try:
    import netCDF4  # noqa: F401
except ImportError:
    import warnings

    warnings.warn("Unable to import netCDF4.", ImportWarning)
    WITH_NETCDF4 = False
else:
    WITH_NETCDF4 = True

__all__ = [
    "read_netcdf",
    "write_netcdf",
    "write_raster_netcdf",
    "NotRasterGridError",
    "WITH_NETCDF4",
]

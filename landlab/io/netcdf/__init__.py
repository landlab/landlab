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

NETCDF4_EXAMPLE_FILE = os.path.join(
    os.path.dirname(__file__), "tests", "data", "test-netcdf4.nc"
)
NETCDF3_64BIT_EXAMPLE_FILE = os.path.join(
    os.path.dirname(__file__), "tests", "data", "test-netcdf3-64bit.nc"
)

__all__ = [
    "read_netcdf",
    "write_netcdf",
    "write_raster_netcdf",
    "NotRasterGridError",
    "WITH_NETCDF4",
    "NETCDF4_EXAMPLE_FILE",
    "NETCDF3_64BIT_EXAMPLE_FILE",
]

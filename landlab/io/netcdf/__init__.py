import os

from .errors import NotRasterGridError
from .read import read_netcdf
from .write import write_netcdf, write_raster_netcdf


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
    "NETCDF4_EXAMPLE_FILE",
    "NETCDF3_64BIT_EXAMPLE_FILE",
]

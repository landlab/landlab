import os

from .read import read_netcdf
from .write import write_netcdf
from .errors import NotRasterGridError

try:
    import netCDF4 as nc4
except ImportError:
    import warnings
    warnings.warn('Unable to import netCDF4.', ImportWarning)
    WITH_NETCDF4 = False
else:
    WITH_NETCDF4 = True

NETCDF4_EXAMPLE_FILE = os.path.join(os.path.dirname(__file__), 'tests', 'data',
                                    'test-netcdf4.nc')
NETCDF3_64BIT_EXAMPLE_FILE = os.path.join(os.path.dirname(__file__), 'tests',
                                          'data', 'test-netcdf3-64bit.nc')

__all__ = ('read_netcdf', 'write_netcdf', 'NotRasterGridError',
           'WITH_NETCDF4', 'NETCDF4_EXAMPLE_FILE',
           'NETCDF3_64BIT_EXAMPLE_FILE')

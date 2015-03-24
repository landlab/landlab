try:
    import netCDF4 as nc4
except ImportError:
    import warnings
    warnings.warn('Unable to import netCDF4.', ImportWarning)
    WITH_NETCDF4 = False
else:
    WITH_NETCDF4 = True


from .read import read_netcdf
from .write import write_netcdf
from .errors import NotRasterGridError


__all__ = ['read_netcdf', 'write_netcdf', 'NotRasterGridError',
           'WITH_NETCDF4', ]

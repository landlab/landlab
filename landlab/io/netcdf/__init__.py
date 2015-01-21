try:
    import netCDF4 as nc4
except ImportError:
    import warnings
    warnings.warn('Unable to import netCDF4.', ImportWarning)
    WITH_NETCDF4 = False
else:
    WITH_NETCDF4 = True


from landlab.io.netcdf.read import read_netcdf
from landlab.io.netcdf.write import write_netcdf
from landlab.io.netcdf.errors import NotRasterGridError

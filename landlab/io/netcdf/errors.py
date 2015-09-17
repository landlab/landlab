#! /usr/bin/env


class Error(Exception):
    """
    Base class for errors in this package.
    """
    pass


class NotRasterGridError(Error):
    """
    Raise this error if the grid defined in the netcdf file is not
    uniform rectilinear with constant spacing in all dimensions.
    """
    pass

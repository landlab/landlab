#! /usr/bin/env
"""Exceptions to raise for the netcdf module."""


class Error(Exception):

    """Base class for errors in this package."""

    pass


class NotRasterGridError(Error):

    """Raise if grid is not uniform rectilinear.

    Raise this error if the grid defined in the netcdf file is not
    uniform rectilinear with constant spacing in all dimensions.
    """

    pass

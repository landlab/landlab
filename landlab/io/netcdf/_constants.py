#!/usr/bin/env python
"""Constants used with the netcdf module."""


_DIMENSION_NAMES = ["ni", "nj", "nk"]
_AXES_NAMES = ["x", "y", "z"]

_NP_TO_NC_TYPE = {
    "float32": "f4",
    "float64": "f8",
    "int8": "i1",
    "int16": "i2",
    "int32": "i4",
    "int64": "i8",
    "uint8": "u1",
    "uint16": "u2",
    "uint32": "u4",
    "uint64": "u8",
    "bool": "i1",
}


_AXIS_DIMENSION_NAMES = ["nk", "nj", "ni"]
_AXIS_COORDINATE_NAMES = ["z", "y", "x"]

_DIMENSION_NAMES = set(_AXIS_DIMENSION_NAMES + ["nt"])
_COORDINATE_NAMES = set(_AXIS_COORDINATE_NAMES + ["t"])


_GRID_MAPPING_VARIABLES = [
    "crs_wkt",
    "earth_radius",
    "false_easting",
    "false_northing",
    "grid_mapping_name",
    "grid_north_pole_latitude",
    "grid_north_pole_longitude",
    "inverse_flattening",
    "latitude_of_projection_origin",
    "longitude_of_central_meridian",
    "longitude_of_prime_meridian",
    "longitude_of_projection_origin",
    "north_pole_grid_longitude",
    "perspective_point_height",
    "scale_factor_at_central_meridian",
    "scale_factor_at_projection_origin",
    "semi_major_axis",
    "semi_minor_axis",
    "standard_parallel",
    "straight_vertical_longitude_from_pole",
]

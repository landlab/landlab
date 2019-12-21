#! /usr/bin/env python
import os

import numpy as np

from landlab import RasterModelGrid
from landlab.io import read_esri_ascii
from landlab.io.netcdf import read_netcdf


def test_save_esri_ascii(tmpdir):
    grid = RasterModelGrid((4, 5), xy_spacing=2.0)
    grid.add_field("node", "air__temperature", np.arange(20.0))

    with tmpdir.as_cwd():
        grid.save("test.asc", format="esri-ascii")
        assert os.path.isfile("test.asc")


def test_add_extension(tmpdir):
    grid = RasterModelGrid((4, 5), xy_spacing=2.0)
    grid.add_field("node", "air__temperature", np.arange(20.0))

    with tmpdir.as_cwd():
        grid.save("test", format="esri-ascii")
        assert os.path.isfile("test.asc")

    with tmpdir.as_cwd():
        grid.save("test", format="netcdf")
        assert os.path.isfile("test.nc")


def test_replace_extension(tmpdir):
    grid = RasterModelGrid((4, 5), xy_spacing=2.0)
    grid.add_field("node", "air__temperature", np.arange(20.0))

    with tmpdir.as_cwd():
        grid.save("test.nc", format="esri-ascii")
        assert os.path.isfile("test.asc")

    with tmpdir.as_cwd():
        grid.save("test.asc", format="netcdf")
        assert os.path.isfile("test.nc")


def test_guess_format(tmpdir):
    grid = RasterModelGrid((4, 5), xy_spacing=2.0)
    grid.add_field("node", "air__temperature", np.arange(20.0))

    with tmpdir.as_cwd():
        grid.save("test.asc")
        assert os.path.isfile("test.asc")
        read_esri_ascii("test.asc")

    with tmpdir.as_cwd():
        grid.save("test.nc")
        assert os.path.isfile("test.nc")
        read_netcdf("test.nc")


def test_names_keyword_as_str(tmpdir):
    grid = RasterModelGrid((4, 5), xy_spacing=2.0)
    grid.add_field("air__temperature", np.arange(20.0), at="node")
    grid.add_field("land_surface__elevation", np.arange(20.0), at="node")

    with tmpdir.as_cwd():
        grid.save("test.asc", names="land_surface__elevation")
        assert os.path.isfile("test.asc")
        read_esri_ascii("test.asc")


def test_names_keyword_as_list(tmpdir):
    grid = RasterModelGrid((4, 5), xy_spacing=2.0)
    grid.add_field("air__temperature", np.arange(20.0), at="node")
    grid.add_field("land_surface__elevation", np.arange(20.0), at="node")

    with tmpdir.as_cwd():
        grid.save("test.asc", names=["land_surface__elevation", "air__temperature"])
        files = ["test_land_surface__elevation.asc", "test_air__temperature.asc"]
        for fname in files:
            assert os.path.isfile(fname)
            read_esri_ascii(fname)

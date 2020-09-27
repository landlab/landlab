#!/usr/bin/env python3

"""unit tests for landlab.io.obj module"""
import pathlib

import pytest

from landlab import HexModelGrid, RasterModelGrid
from landlab.io import write_obj

LITTLE_HEX_OBJ = """# landlabgrid
#
g landlabgrid
v 1.0 0.0 0.0
v 3.0 0.0 0.0
v 0.0 1.732051 0.0
v 2.0 1.732051 1.0
v 4.0 1.732051 0.0
v 1.0 3.464102 0.0
v 3.0 3.464102 0.0
f 4// 1// 2//
f 4// 3// 1//
f 5// 4// 2//
f 6// 3// 4//
f 7// 4// 5//
f 7// 6// 4//
"""

LITTLE_RAST_OBJ = """# landlabgrid
#
g landlabgrid
v 0.0 0.0 0.0
v 2.0 0.0 0.0
v 4.0 0.0 0.0
v 0.0 2.0 0.0
v 2.0 2.0 1.0
v 4.0 2.0 0.0
v 0.0 4.0 0.0
v 2.0 4.0 0.0
v 4.0 4.0 0.0
f 5// 4// 1//
f 1// 2// 5//
f 6// 5// 2//
f 2// 3// 6//
f 8// 7// 4//
f 4// 5// 8//
f 9// 8// 5//
f 5// 6// 9//
"""


def test_write_to_filelike(tmpdir):
    grid = RasterModelGrid((3, 3), xy_spacing=2.0)
    z = grid.add_zeros("topographic__elevation", at="node")
    z[4] = 1.0

    with tmpdir.as_cwd():
        with open("test_quad.obj", "w") as fp:
            write_obj(fp, grid)
        with open("test_quad.obj", "r") as fp:
            assert fp.read() == LITTLE_RAST_OBJ


@pytest.mark.parametrize("fname", (pathlib.Path("test_hex.obj"), "test_hex.obj"))
def test_write_hex_to_path(tmpdir, fname):
    grid = HexModelGrid((3, 2), spacing=2.0)
    z = grid.add_zeros("topographic__elevation", at="node")
    z[3] = 1.0

    with tmpdir.as_cwd():
        actual = write_obj(fname, grid)
        assert actual == fname
        with open(fname, "r") as fp:
            assert fp.read() == LITTLE_HEX_OBJ


def test_write_raster(tmpdir):
    grid = RasterModelGrid((3, 3), xy_spacing=2.0)
    z = grid.add_zeros("topographic__elevation", at="node")
    z[4] = 1.0

    with tmpdir.as_cwd():
        write_obj("test_quad.obj", grid)
        with open("test_quad.obj", "r") as fp:
            assert fp.read() == LITTLE_RAST_OBJ


def test_field_name(tmpdir):
    grid = RasterModelGrid((3, 3), xy_spacing=2.0)
    z = grid.add_zeros("z", at="node")
    z[4] = 1.0

    with tmpdir.as_cwd():
        write_obj("test_quad.obj", grid, field_for_z="z")
        with open("test_quad.obj", "r") as fp:
            assert fp.read() == LITTLE_RAST_OBJ


def test_clobber(tmpdir):
    grid = RasterModelGrid((3, 3), xy_spacing=2.0)
    z = grid.add_zeros("topographic__elevation", at="node")
    z[4] = 1.0

    with tmpdir.as_cwd():
        with open("test_quad.obj", "w") as fp:
            pass

        with pytest.raises(ValueError):
            write_obj("test_quad.obj", grid)

        write_obj("test_quad.obj", grid, clobber=True)
        with open("test_quad.obj", "r") as fp:
            assert fp.read() == LITTLE_RAST_OBJ

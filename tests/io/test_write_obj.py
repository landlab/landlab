#!/usr/bin/env python3

"""unit tests for landlab.io.obj module"""

import os
import tempfile
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


def remove_newline_chars(mystr):
    mystr = mystr.replace("\n", "")
    mystr = mystr.replace("\r", "")
    return mystr


def test_write_obj():

    grid = HexModelGrid((3, 2), spacing=2.0)
    z = grid.add_zeros("topographic__elevation", at="node")
    z[3] = 1.0
    with tempfile.TemporaryDirectory() as tmpdirname:
        fname = os.path.join(tmpdirname, "test_tri.obj")
        write_obj(fname, grid)
        file = open(fname)
        contents = file.read()
        contents = remove_newline_chars(contents)
        assert contents == remove_newline_chars(LITTLE_HEX_OBJ)

    grid = RasterModelGrid((3, 3), xy_spacing=2.0)
    z = grid.add_zeros("topographic__elevation", at="node")
    z[4] = 1.0
    with tempfile.TemporaryDirectory() as tmpdirname:
        fname = os.path.join(tmpdirname, "test_quad.obj")
        write_obj(fname, grid)
        file = open(fname)
        contents = file.read()
        contents = remove_newline_chars(contents)
        assert contents == remove_newline_chars(LITTLE_RAST_OBJ)

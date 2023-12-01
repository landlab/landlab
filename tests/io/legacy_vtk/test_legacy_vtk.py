#! /usr/bin/env python

import io

import numpy as np

from landlab import RasterModelGrid
from landlab.io import write_legacy_vtk

EXPECTED_VTK_FOR_RASTER = """# vtk DataFile Version 2.0
Landlab output
ASCII
DATASET UNSTRUCTURED_GRID

POINTS 9 float
0.0 0.0 0.0
1.0 0.0 0.0
2.0 0.0 0.0
0.0 1.0 0.0
1.0 1.0 1.0
2.0 1.0 0.0
0.0 2.0 0.0
1.0 2.0 0.0
2.0 2.0 0.0

CELLS 4 20
4 4 3 0 1
4 5 4 1 2
4 7 6 3 4
4 8 7 4 5

CELL_TYPES 4
9
9
9
9

POINT_DATA 9
SCALARS topographic__elevation float 1
LOOKUP_TABLE default
0.0
0.0
0.0
0.0
1.0
0.0
0.0
0.0
0.0
SCALARS surface_water__depth float 1
LOOKUP_TABLE default
0.0
1.0
2.0
3.0
4.0
5.0
6.0
7.0
8.0

"""


def test_raster_grid():
    # Create a tiny grid with 1 core node and 6 boundary nodes
    grid = RasterModelGrid((3, 3))

    # Add two fields with made-up values
    topo = grid.add_zeros("topographic__elevation", at="node")
    topo[4] = 1.0
    water = grid.add_zeros("surface_water__depth", at="node")
    water[:] = np.arange(len(water))

    # Write output in legacy VTK format
    vtk_file = write_legacy_vtk(io.StringIO(), grid)

    assert vtk_file.getvalue() == EXPECTED_VTK_FOR_RASTER

import io

import numpy as np

import landlab.io.legacy_vtk as vtk
from landlab import RasterModelGrid

EXPECTED_VTK_FOR_RASTER = """\
# vtk DataFile Version 2.0
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
0.0"""


def test_raster_grid():
    # Create a tiny grid with 1 core node and 6 boundary nodes
    grid = RasterModelGrid((3, 3))

    # Add two fields with made-up values
    topo = grid.add_zeros("topographic__elevation", at="node")
    topo[4] = 1.0
    water = grid.add_zeros("surface_water__depth", at="node")
    water[:] = np.arange(len(water))

    # Write output in legacy VTK format
    vtk_file = vtk.write_legacy_vtk(io.StringIO(), grid)

    assert vtk_file.getvalue() == EXPECTED_VTK_FOR_RASTER


def test_dump(tmpdir):
    grid = RasterModelGrid((3, 3))

    topo = grid.add_zeros("topographic__elevation", at="node")
    topo[4] = 1.0
    water = grid.add_zeros("surface_water__depth", at="node")
    water[:] = np.arange(len(water))

    actual = vtk.dump(grid, z_coord=topo)

    assert actual == EXPECTED_VTK_FOR_RASTER

    with tmpdir.as_cwd():
        with open("grid.vtk", "w") as fp:
            vtk.dump(grid, z_coord=topo, stream=fp)
        with open("grid.vtk") as fp:
            actual = fp.read()

    assert actual == EXPECTED_VTK_FOR_RASTER


def test_dump_patches():
    grid = RasterModelGrid((3, 4))

    grid.at_patch["topographic__elevation"] = np.arange(grid.number_of_patches)
    grid.at_patch["surface_water__depth"] = np.arange(grid.number_of_patches)

    actual = vtk.dump(grid)

    assert ("CELL_DATA 6" in actual) and ("POINT_DATA" not in actual)


def test_dump_include():
    grid = RasterModelGrid((3, 4))

    grid.at_patch["topographic__elevation"] = np.arange(grid.number_of_patches)
    grid.at_patch["surface_water__depth"] = np.arange(grid.number_of_patches)

    actual = vtk.dump(grid, include="*elevation*")

    assert ("CELL_DATA 6" in actual) and ("POINT_DATA" not in actual)
    assert ("topographic__elevation" in actual) and (
        "surface_water__depth" not in actual
    )


def test_dump_exclude():
    grid = RasterModelGrid((3, 4))

    grid.at_patch["topographic__elevation"] = np.arange(grid.number_of_patches)
    grid.at_patch["surface_water__depth"] = np.arange(grid.number_of_patches)

    actual = vtk.dump(grid, exclude="*")
    assert ("CELL_DATA" not in actual) and ("POINT_DATA" not in actual)

    actual = vtk.dump(grid, exclude="*water*")

    assert f"CELL_DATA {grid.number_of_patches}" in actual
    assert "POINT_DATA" not in actual
    assert ("topographic__elevation" in actual) and (
        "surface_water__depth" not in actual
    )


def test_dump_dual():
    grid = RasterModelGrid((3, 4))

    grid.at_cell["topographic__elevation"] = np.arange(grid.number_of_cells)
    grid.at_cell["surface_water__depth"] = np.arange(grid.number_of_cells)

    grid.at_corner["topographic__elevation"] = np.arange(grid.number_of_corners)
    grid.at_corner["surface_water__depth"] = np.arange(grid.number_of_corners)

    actual = vtk.dump(grid, at="corner")

    assert f"POINTS {grid.number_of_corners}" in actual
    assert f"POINT_DATA {grid.number_of_corners}" in actual
    assert f"CELLS {grid.number_of_cells}" in actual
    assert f"CELL_DATA {grid.number_of_cells}" in actual

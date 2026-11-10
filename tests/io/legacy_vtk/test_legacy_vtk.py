import io

import numpy as np
import pytest
from numpy.testing import assert_array_equal

import landlab.io.legacy_vtk as vtk
from landlab import RasterModelGrid
from landlab.io.legacy_vtk import _format_vtk_scalar_data

EXPECTED_VTK_FOR_RASTER = """\
# vtk DataFile Version 2.0
Landlab output
ASCII
DATASET UNSTRUCTURED_GRID

POINTS 9 double
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

SCALARS surface_water__depth double 1
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

SCALARS topographic__elevation double 1
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


def test_write_tiny_node_field():
    grid = RasterModelGrid((3, 3))
    grid.add_zeros("topographic__elevation", at="node")
    deposition = grid.add_zeros("landslide__deposition", at="node")
    deposition[:5] = [
        1.4214660673135646e-38,
        2.0064221222504523e-39,
        3.6073920684198113e-40,
        5.434904772579652e-41,
        4.089438831653804e-42,
    ]

    output = vtk.write_legacy_vtk(io.StringIO(), grid).getvalue()
    scalar_data = output.split("SCALARS landslide__deposition double 1\n")[1]
    values = scalar_data.splitlines()[1 : grid.number_of_nodes + 1]

    assert_array_equal(np.asarray(values, dtype=float), deposition)


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


TINY = np.finfo(np.double).tiny


@pytest.mark.parametrize(
    "values, expected",
    (
        ([-TINY * 0.5, -TINY, -TINY * 2.0], [0.0, -TINY, -TINY * 2.0]),
        ([TINY * 0.5, TINY, TINY * 2.0], [0.0, TINY, TINY * 2.0]),
        ([-1.0, 0.0, 1.0], [-1.0, 0.0, 1.0]),
        ([-1, 0, 1], [-1.0, 0.0, 1.0]),
    ),
)
def test_format_vtk_scalar_data_with_tiny_data(values, expected):
    actual = _format_vtk_scalar_data(values, name="foobar").splitlines()

    assert actual[0] == "SCALARS foobar double 1"
    assert actual[1] == "LOOKUP_TABLE default"
    assert actual[2:] == list(map(str, expected))
    assert_array_equal(np.asarray(list(map(float, actual[2:]))), expected)


def test_format_vtk_scalar_data_values_not_modified():
    initial = np.array([-1.0, -TINY * 0.5, 0.0, TINY * 0.5, 1.0])
    values = np.copy(initial)

    _format_vtk_scalar_data(values, name="foobar")

    assert_array_equal(values, initial)

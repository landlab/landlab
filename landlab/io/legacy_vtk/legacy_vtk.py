#! /usr/bin/env python

import os
import pathlib


def _write_vtk_header(file_like):
    """Write the file header."""
    file_like.write("# vtk DataFile Version 2.0\n")
    file_like.write("Landlab output\n")
    file_like.write("ASCII\n")
    file_like.write("DATASET UNSTRUCTURED_GRID\n\n")


def _write_vtk_points(grid, file_like, z_at_node):
    """Write the POINTS section of the file(-like)"""
    x = grid.x_of_node
    y = grid.y_of_node
    z = z_at_node
    file_like.write("POINTS " + str(grid.number_of_nodes) + " float\n")
    for i in range(grid.number_of_nodes):
        file_like.write(str(x[i]) + " " + str(y[i]) + " " + str(z[i]) + "\n")
    file_like.write("\n")


def _write_vtk_patches(grid, file_like):
    """Write the CELLS section (in a Landlab grid these are patches)"""
    num_patches = grid.number_of_patches
    nodes_per_patch = len(grid.nodes_at_patch[0])
    file_like.write(
        "CELLS "
        + str(num_patches)
        + " "
        + str((nodes_per_patch + 1) * num_patches)
        + "\n"
    )
    for i in range(grid.number_of_patches):
        file_like.write(str(nodes_per_patch))
        for j in range(nodes_per_patch):
            file_like.write(" " + str(grid.nodes_at_patch[i, j]))
        file_like.write("\n")
    file_like.write("\n")


def _write_vtk_cell_types(grid, file_like):
    """Write the CELL_TYPES section (triangles or quads)"""
    file_like.write("CELL_TYPES " + str(grid.number_of_patches) + "\n")
    if len(grid.nodes_at_patch[0]) == 3:  # triangles
        cell_type = "5\n"  # vtk code for a triangle
    else:
        cell_type = "9\n"  # vtk code for a quad
    for _ in range(grid.number_of_patches):
        file_like.write(cell_type)
    file_like.write("\n")


def _write_scalar_data(grid, file_like, field):
    """Write the SCALARS section for a given field"""
    file_like.write("SCALARS " + field + " float 1\n")
    file_like.write("LOOKUP_TABLE default\n")
    for i in range(grid.number_of_nodes):
        file_like.write(str(grid.at_node[field][i]))
        file_like.write("\n")


def _write_vtk_point_data(grid, file_like, fields):
    """Write the POINT_DATA section, which in turn writes a SCALARS
    section for each field in `fields`"""
    file_like.write("POINT_DATA " + str(grid.number_of_nodes) + "\n")
    for fieldname in fields:
        _write_scalar_data(grid, file_like, fieldname)
    file_like.write("\n")


def write_legacy_vtk(
    path, grid, z_at_node="topographic__elevation", fields=None, clobber=False
):
    """
    Write grid and field to a legacy VTK format file or file-like object.

    Parameters
    ----------
    path : file-like
        Path to file or a file-like object
    grid : Landlab grid object
        The grid for which to output data
    z_at_node : str or (n_nodes, ) ndarray
        Field name or array of values to use for z coordinate
    fields : list of str (optional)
        List of node fields to output; default is all node fields
    clobber : bool (optional)
        Ok to overwrite existing file (default False)

    Examples
    --------
    >>> import io
    >>> import numpy as np
    >>> from landlab import HexModelGrid
    >>> from landlab.io.legacy_vtk import write_legacy_vtk

    >>> grid = HexModelGrid((3, 2))
    >>> topo = grid.add_zeros("topographic__elevation", at="node")
    >>> topo[:] = np.arange(len(topo))
    >>> water = grid.add_zeros("surface_water__depth", at="node")
    >>> water[:] = 0.1 * (7.0 - topo)

    >>> vtk_file = write_legacy_vtk(io.StringIO(), grid)
    >>> lines = vtk_file.getvalue().splitlines()
    >>> print(lines[0])
    # vtk DataFile Version 2.0
    >>> for i in range(5, 13):
    ...     print(lines[i])
    POINTS 7 float
    0.5 0.0 0.0
    1.5 0.0 1.0
    0.0 0.866025 2.0
    1.0 0.866025 3.0
    2.0 0.866025 4.0
    0.5 1.732051 5.0
    1.5 1.732051 6.0
    >>> for i in range(14, 21):
    ...     print(lines[i])
    CELLS 6 24
    3 3 0 1
    3 3 2 0
    3 4 3 1
    3 5 2 3
    3 6 3 4
    3 6 5 3
    >>> for i in range(22, 29):
    ...     print(lines[i])
    CELL_TYPES 6
    5
    5
    5
    5
    5
    5
    >>> for i in range(30, 49):
    ...     print(lines[i])
    POINT_DATA 7
    SCALARS topographic__elevation float 1
    LOOKUP_TABLE default
    0.0
    1.0
    2.0
    3.0
    4.0
    5.0
    6.0
    SCALARS surface_water__depth float 1
    LOOKUP_TABLE default
    0.7
    0.6
    0.5
    0.4
    0.3
    0.2
    0.1
    """
    if isinstance(z_at_node, str):
        z_at_node = grid.at_node[z_at_node]
    if fields is None:
        fields = grid.at_node.keys()

    if isinstance(path, (str, pathlib.Path)):
        if os.path.exists(path) and not clobber:
            raise ValueError(f"file exists ({path})")

        with open(path, "w") as fp:
            _write_legacy_vtk_to_filelike(fp, grid, z_at_node, fields)
    else:
        _write_legacy_vtk_to_filelike(path, grid, z_at_node, fields)

    return path


def _write_legacy_vtk_to_filelike(file_like, grid, z_at_node, fields):
    """Write output to specified file_like"""
    _write_vtk_header(file_like)
    _write_vtk_points(grid, file_like, z_at_node)
    _write_vtk_patches(grid, file_like)
    _write_vtk_cell_types(grid, file_like)
    _write_vtk_point_data(grid, file_like, fields)

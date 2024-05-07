import os
import pathlib
from enum import IntEnum
from enum import unique

import numpy as np


def dump(grid, stream=None, include="*", exclude=None, z_coord=0.0, at="node"):
    """Format a grid in VTK legacy format.

    Parameters
    ----------
    grid : ModelGrid
        A *Landlab* grid.
    stream : file_like, optional
        File-like object to write the formatted output. If not provided,
        return the formatted vtk as a string.
    include : str, or iterable of str, optional
        Glob-style pattern for field names to include.
    exclude : str, or iterable of str, optional
        Glob-style pattern for field names to exclude.
    z_coord : array_like or str, optional
        If the grid does not have a *z* coordinate, use this value. If
        a ``str``, use the corresponding field.
    at : 'node' or 'corner', optional
        Use either the grid's *node*s (and *patches*) or *corners* (and *cells*).

    Returns
    -------
    str or ``None``
        The grid formatted as legacy VTK or ``None`` if an output stream
        was provided.

    Examples
    --------
    >>> import os
    >>> import numpy as np
    >>> from landlab import HexModelGrid
    >>> import landlab.io.legacy_vtk.legacy_vtk as vtk

    >>> grid = HexModelGrid((3, 2))
    >>> topo = grid.add_zeros("topographic__elevation", at="node")
    >>> topo[:] = np.arange(len(topo))
    >>> water = grid.add_zeros("surface_water__depth", at="node")
    >>> water[:] = (7.0 - topo) / 10.0

    >>> lines = vtk.dump(grid).splitlines()
    >>> print(os.linesep.join(lines[:4]))
    # vtk DataFile Version 2.0
    Landlab output
    ASCII
    DATASET UNSTRUCTURED_GRID
    """
    content = _format_as_vtk(
        grid, include=include, exclude=exclude, z_coord=z_coord, at=at
    )
    if stream is None:
        return content
    else:
        stream.write(content)


def _format_as_vtk(grid, include="*", exclude=None, z_coord=0.0, at="node"):
    if at not in ("node", "corner"):
        raise ValueError(f"`at` keyword must be one of 'node' or 'corner' ({at})")

    if at == "node":
        point, cell = "node", "patch"
    else:
        point, cell = "corner", "cell"

    at_point = getattr(grid, f"at_{point}")
    at_cell = getattr(grid, f"at_{cell}")
    points_at_cell = getattr(grid, f"{point}s_at_{cell}")

    if isinstance(z_coord, str):
        z_coord = at_point[z_coord]

    try:
        coords_of_point = getattr(grid, f"coords_of_{point}")
    except AttributeError:
        coords_of_point = np.pad(getattr(grid, f"xy_of_{point}"), ((0, 0), (0, 1)))
        coords_of_point[:, -1] = z_coord

    content = [
        _format_vtk_header(),
        _format_vtk_points(coords_of_point),
        _format_vtk_cells(points_at_cell),
    ]

    fields = grid.fields(include=include, exclude=exclude)

    point_fields = {
        field.split(":")[1]: at_point[field.split(":")[1]]
        for field in fields if field.startswith(f"at_{point}:")
    }
    if point_fields:
        content.append(_format_vtk_point_data(point_fields))

    cell_fields = {
        field.split(":")[1]: at_point[field.split(":")[1]]
        for field in fields if field.startswith(f"at_{cell}:")
    }
    if cell_fields:
        content.append(_format_vtk_cell_data(cell_fields))

    return (2 * os.linesep).join(content)


@unique
class CellType(IntEnum):
    TRIANGLE = 3
    QUAD = 9
    POLYGON = 7


VTK_CELL_TYPE = {
    3: CellType.TRIANGLE,
    4: CellType.QUAD,
}


def _format_vtk_header():
    return """\
# vtk DataFile Version 2.0
Landlab output
ASCII
DATASET UNSTRUCTURED_GRID"""


def _format_vtk_points(coords_of_points):
    return os.linesep.join(
        [f"POINTS {len(coords_of_points)} float"]
        + [
            " ".join(str(coord) for coord in coords)
            for coords in coords_of_points
        ]
    )


def _format_vtk_cells(points_at_cell):
    cells = []
    types = []
    n_points = 0
    for cell in points_at_cell:
        points = [point for point in cell if point >= -1]
        if points:
            cells.append(f"{len(points)} " + " ".join(str(point) for point in points))
            types.append(str(VTK_CELL_TYPE.get(len(points), CellType.POLYGON)))
            n_points += len(points) + 1

    cells_section = os.linesep.join(
        [f"CELLS {len(cells)} {n_points}"] + cells
    )

    types_section = os.linesep.join(
        [f"CELL_TYPES {len(types)}"] + types
    )

    return (2 * os.linesep).join([cells_section, types_section])


def _format_vtk_point_data(point_data):
    content = []
    for name, value_at_point in point_data.items():
        content.append(_format_vtk_scalar_data(value_at_point, name=name))
        number_of_points = len(value_at_point)

    return (2 * os.linesep).join([f"POINT_DATA {number_of_points}"] + content)


def _format_vtk_cell_data(cell_data):
    content = []
    for name, value_at_cell in cell_data.items():
        content.append(_format_vtk_scalar_data(value_at_cell, name=name))
        number_of_cells = len(value_at_cell)

    return (2 * os.linesep).join([f"CELL_DATA {number_of_cells}"] + content)


def _format_vtk_scalar_data(values, name="data"):
    return os.linesep.join(
        [
            f"""\
SCALARS {name} float 1
LOOKUP_TABLE default"""
        ] + [str(value) for value in values]
    )


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
    file_like.write(f"POINTS {grid.number_of_nodes} float\n")
    for i in range(grid.number_of_nodes):
        file_like.write(f"{x[i]} {y[i]} {z[i]}\n")
    file_like.write("\n")


def _write_vtk_patches(grid, file_like):
    """Write the CELLS section (in a Landlab grid these are patches)"""
    num_patches = grid.number_of_patches
    nodes_per_patch = len(grid.nodes_at_patch[0])
    file_like.write(f"CELLS {num_patches} {(nodes_per_patch + 1) * num_patches}\n")
    for i in range(grid.number_of_patches):
        file_like.write(str(nodes_per_patch))
        for j in range(nodes_per_patch):
            file_like.write(f" {grid.nodes_at_patch[i, j]}")
        file_like.write("\n")
    file_like.write("\n")


def _write_vtk_cell_types(grid, file_like):
    """Write the CELL_TYPES section (triangles or quads)"""
    file_like.write(f"CELL_TYPES {grid.number_of_patches}\n")
    if len(grid.nodes_at_patch[0]) == 3:  # triangles
        cell_type = 5  # vtk code for a triangle
    else:
        cell_type = 9  # vtk code for a quad
    for _ in range(grid.number_of_patches):
        file_like.write(f"{cell_type}\n")

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
    >>> water[:] = (7.0 - topo) / 10.0

    >>> vtk_file = write_legacy_vtk(io.StringIO(), grid)
    >>> lines = vtk_file.getvalue().splitlines()
    >>> print(lines[0])
    # vtk DataFile Version 2.0
    >>> for i in range(5, 13):
    ...     print(lines[i])
    ...
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
    ...
    CELLS 6 24
    3 3 0 1
    3 3 2 0
    3 4 3 1
    3 5 2 3
    3 6 3 4
    3 6 5 3
    >>> for i in range(22, 29):
    ...     print(lines[i])
    ...
    CELL_TYPES 6
    5
    5
    5
    5
    5
    5
    >>> for i in range(30, 49):
    ...     print(lines[i])
    ...
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

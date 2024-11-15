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
    >>> import landlab.io.legacy_vtk as vtk

    >>> grid = HexModelGrid((3, 2))

    >>> topo = np.arange(grid.number_of_nodes)
    >>> grid.at_node["topographic__elevation"] = topo
    >>> grid.at_node["surface_water__depth"] = (7.0 - topo) / 10.0

    >>> lines = vtk.dump(grid, z_coord=topo).splitlines()
    >>> print(os.linesep.join(lines[:4]))
    # vtk DataFile Version 2.0
    Landlab output
    ASCII
    DATASET UNSTRUCTURED_GRID

    The x, y, and z coordinates of each grid node (VTK calls these
    "points")

    >>> print(os.linesep.join(lines[5:13]))
    POINTS 7 float
    0.5 0.0 0.0
    1.5 0.0 1.0
    0.0 0.866025 2.0
    1.0 0.866025 3.0
    2.0 0.866025 4.0
    0.5 1.732051 5.0
    1.5 1.732051 6.0

    Grid nodes that form each patch (VTK calls these "cells").

    >>> print(os.linesep.join(lines[14:21]))
    CELLS 6 24
    3 3 0 1
    3 3 2 0
    3 4 3 1
    3 5 2 3
    3 6 3 4
    3 6 5 3

    The type of each patch. A value of 5 is VTK code for triangle.

    >>> print(os.linesep.join(lines[22:29]))
    CELL_TYPES 6
    5
    5
    5
    5
    5
    5

    The data fields at grid nodes.

    >>> print(os.linesep.join(lines[30:51]))
    POINT_DATA 7
    <BLANKLINE>
    SCALARS surface_water__depth float 1
    LOOKUP_TABLE default
    0.7
    0.6
    0.5
    0.4
    0.3
    0.2
    0.1
    <BLANKLINE>
    SCALARS topographic__elevation float 1
    LOOKUP_TABLE default
    0.0
    1.0
    2.0
    3.0
    4.0
    5.0
    6.0

    To write the dual grid (i.e. corners and cells) use the ``at``
    keyword.

    >>> lines = vtk.dump(grid, at="corner").splitlines()
    >>> print(os.linesep.join(lines[5:12]))
    POINTS 6 float
    1.0 0.288675 0.0
    0.5 0.57735 0.0
    1.5 0.57735 0.0
    0.5 1.154701 0.0
    1.5 1.154701 0.0
    1.0 1.443376 0.0
    """
    content = _format_as_vtk(
        grid, include=include, exclude=exclude, z_coord=z_coord, at=at
    )
    if stream is None:
        return content
    else:
        stream.write(content)


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
    >>> import os
    >>> import numpy as np
    >>> from landlab import HexModelGrid
    >>> from landlab.io.legacy_vtk import write_legacy_vtk

    >>> grid = HexModelGrid((3, 2))

    >>> topo = np.arange(grid.number_of_nodes)
    >>> grid.at_node["topographic__elevation"] = topo
    >>> grid.at_node["surface_water__depth"] = (7.0 - topo) / 10.0

    >>> vtk_file = write_legacy_vtk(io.StringIO(), grid)
    >>> lines = vtk_file.getvalue().splitlines()
    >>> print(lines[0])
    # vtk DataFile Version 2.0

    >>> print(os.linesep.join(lines[5:13]))
    POINTS 7 float
    0.5 0.0 0.0
    1.5 0.0 1.0
    0.0 0.866025 2.0
    1.0 0.866025 3.0
    2.0 0.866025 4.0
    0.5 1.732051 5.0
    1.5 1.732051 6.0

    >>> print(os.linesep.join(lines[14:21]))
    CELLS 6 24
    3 3 0 1
    3 3 2 0
    3 4 3 1
    3 5 2 3
    3 6 3 4
    3 6 5 3

    >>> print(os.linesep.join(lines[22:29]))
    CELL_TYPES 6
    5
    5
    5
    5
    5
    5

    >>> print(os.linesep.join(lines[30:51]))
    POINT_DATA 7
    <BLANKLINE>
    SCALARS surface_water__depth float 1
    LOOKUP_TABLE default
    0.7
    0.6
    0.5
    0.4
    0.3
    0.2
    0.1
    <BLANKLINE>
    SCALARS topographic__elevation float 1
    LOOKUP_TABLE default
    0.0
    1.0
    2.0
    3.0
    4.0
    5.0
    6.0
    """
    if isinstance(z_at_node, str):
        z_at_node = grid.at_node[z_at_node]

    if fields is None:
        fields = grid.at_node.keys()

    if isinstance(fields, str):
        fields = [fields]
    fields = [f"at_node:{field}" for field in fields]

    if isinstance(path, (str, pathlib.Path)):
        if os.path.exists(path) and not clobber:
            raise ValueError(f"file exists ({path})")

        with open(path, "w") as fp:
            dump(grid, stream=fp, include=fields, z_coord=z_at_node, at="node")

    else:
        dump(grid, stream=path, include=fields, z_coord=z_at_node, at="node")

    return path


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
        coords_of_point = getattr(grid, f"xy_of_{point}")

    if coords_of_point.shape[1] == 2:
        coords_of_point = np.pad(getattr(grid, f"xy_of_{point}"), ((0, 0), (0, 1)))
        try:
            coords_of_point[:, -1] = z_coord
        except ValueError as e:
            e.add_note(
                f"The grid has {len(coords_of_point)} {at}s but the provided value"
                f" value for z_coord has size {np.shape(np.atleast_1d(z_coord))}."
            )
            raise

    content = [
        _format_vtk_header(),
        _format_vtk_points(coords_of_point),
        _format_vtk_cells(points_at_cell),
    ]

    fields = grid.fields(include=include, exclude=exclude)

    point_fields = {
        field.split(":")[1]: at_point[field.split(":")[1]]
        for field in sorted(fields)
        if field.startswith(f"at_{point}:")
    }
    if point_fields:
        content.append(_format_vtk_point_data(point_fields))

    cell_fields = {
        field.split(":")[1]: at_cell[field.split(":")[1]]
        for field in sorted(fields)
        if field.startswith(f"at_{cell}:")
    }
    if cell_fields:
        content.append(_format_vtk_cell_data(cell_fields))

    return (2 * "\n").join(content)


@unique
class CellType(IntEnum):
    TRIANGLE = 5
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
    return "\n".join(
        [f"POINTS {len(coords_of_points)} float"]
        + [" ".join(str(coord) for coord in coords) for coords in coords_of_points]
    )


def _format_vtk_cells(points_at_cell):
    cells = []
    types = []
    n_points = 0
    for cell in points_at_cell:
        points = [point for point in cell if point >= -1]
        if points:
            cells.append(f"{len(points)} " + " ".join(str(point) for point in points))
            types.append(format(VTK_CELL_TYPE.get(len(points), CellType.POLYGON)))
            n_points += len(points) + 1

    cells_section = "\n".join([f"CELLS {len(cells)} {n_points}"] + cells)

    types_section = "\n".join([f"CELL_TYPES {len(types)}"] + types)

    return (2 * "\n").join([cells_section, types_section])


def _format_vtk_point_data(point_data):
    content = []
    for name, value_at_point in point_data.items():
        content.append(_format_vtk_scalar_data(value_at_point, name=name))
        number_of_points = len(value_at_point)

    return (2 * "\n").join([f"POINT_DATA {number_of_points}"] + content)


def _format_vtk_cell_data(cell_data):
    content = []
    for name, value_at_cell in cell_data.items():
        content.append(_format_vtk_scalar_data(value_at_cell, name=name))
        number_of_cells = len(value_at_cell)

    return (2 * "\n").join([f"CELL_DATA {number_of_cells}"] + content)


def _format_vtk_scalar_data(values, name="data"):
    return "\n".join(
        [
            f"""\
SCALARS {name} float 1
LOOKUP_TABLE default"""
        ]
        + [str(float(value)) for value in values]
    )

#! /usr/bin/env python
"""Write (x,y,z) data from a Landlab grid + 1 field to a Wavefront OBJ file.

OBJ functions
+++++++++++++

.. autosummary::

    ~landlab.io.obj.dump
"""
from __future__ import annotations

import os
import pathlib
from collections.abc import Iterable
from typing import TextIO

import numpy as np
from numpy.typing import ArrayLike

from landlab.grid.base import ModelGrid


def dump(
    grid: ModelGrid,
    stream: TextIO | None = None,
    values: str | ArrayLike = 0.0,
    at: str = "node",
) -> str | None:
    """Format a grid in Wavefront OBJ format.

    Write grid data Wavefront OBJ format [1]_.

    Parameters
    ----------
    grid : ModelGrid
        A *Landlab* grid.
    stream : file_like, optional
        File-like object to write the formatted output. If not provided,
        return the formatted grid as a string.
    values : array_like or str, optional
        If the grid does not have a *z* coordinate, use this value. If
        a ``str``, use the corresponding field.
    at : 'node' or 'corner', optional
        Use either the grid's *node*s (and *patches*) or *corners* (and *cells*).

    Returns
    -------
    str or ``None``
        The grid formatted as OBJ or ``None`` if an output stream
        was provided.

    References
    ----------

    .. [1] https://www.loc.gov/preservation/digital/formats/fdd/fdd000507.shtml

    Examples
    --------
    >>> from landlab import HexModelGrid
    >>> from landlab.io import obj

    >>> grid = HexModelGrid((3, 2), spacing=2.0)

    Define values at both *nodes* and *corners*.

    >>> z_at_node = grid.add_zeros("topographic__elevation", at="node")
    >>> z_at_node[3] = 1.0

    >>> z_at_corner = grid.add_zeros("topographic__elevation", at="corner")
    >>> z_at_corner[3] = 10.0

    Format using the grid's *nodes*/*patches*.

    >>> print(obj.dump(grid, at="node", values="topographic__elevation"))
    # landlabgrid
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

    Format using the grid's *nodes*/*cells*.

    >>> print(obj.dump(grid, at="corner", values="topographic__elevation"))
    # landlabgrid
    #
    g landlabgrid
    v 2.0 0.57735 0.0
    v 1.0 1.154701 0.0
    v 3.0 1.154701 0.0
    v 1.0 2.309401 10.0
    v 3.0 2.309401 0.0
    v 2.0 2.886751 0.0
    f 5// 6// 4// 2// 1// 3//
    """
    if at not in ("node", "corner"):
        raise ValueError(f"`at` keyword must be one of 'node' or 'corner' ({at!r})")

    if isinstance(values, str):
        z_at_vertex = getattr(grid, f"at_{at}")[values]
    else:
        z_at_vertex = values

    content = _format_as_obj(grid, z_coord=z_at_vertex, at=at)

    if stream is None:
        return content
    else:
        stream.write(content)
        return None


def _format_as_obj(
    grid: ModelGrid,
    z_coord: ArrayLike = 0.0,
    at: str = "node",
) -> str:
    if at == "node":
        vertex, face = "node", "patch"
    else:
        vertex, face = "corner", "cell"

    try:
        coords_of_vertex = getattr(grid, f"coords_of_{vertex}")
    except AttributeError:
        coords_of_vertex = getattr(grid, f"xy_of_{vertex}")

    if coords_of_vertex.shape[1] == 2:
        coords_of_vertex = np.pad(getattr(grid, f"xy_of_{vertex}"), ((0, 0), (0, 1)))
        try:
            coords_of_vertex[:, -1] = z_coord
        except ValueError as e:
            e.add_note(
                f"The grid has {len(coords_of_vertex)} {at}s but the provided value"
                f" for `z_coord` has size {np.shape(np.atleast_1d(z_coord))}."
            )
            raise

    header = """\
# landlabgrid
#
g landlabgrid"""

    return (2 * "\n").join(
        [
            header,
            _format_vertices(coords_of_vertex),
            _format_faces(getattr(grid, f"{vertex}s_at_{face}")),
        ]
    )


def _format_vertices(coords_of_vertices: Iterable[tuple[float, float, float]]) -> str:
    return "\n".join([f"v {x} {y} {z}" for x, y, z in coords_of_vertices])


def _format_faces(vertices_at_face: Iterable[tuple[int, int, int]]) -> str:
    return "\n".join(
        [
            "f " + " ".join(f"{vertex + 1}//" for vertex in vertices if vertex >= 0)
            for vertices in vertices_at_face
        ]
    )


def _write_obj_vertices(file_like, xyz_at_vertex):
    """Write node (x,y,z) coordinates as OBJ vertices.

    Parameters
    ----------
    file_like : file-like
        Opened file-like object to write vertices to.
    xyz_at_vertex : array-like, shape *(n_vertices, 3)*
        (x, y, z) values for each vertex.
    """
    for x, y, z in xyz_at_vertex:
        print(f"v {x} {y} {z}", file=file_like)


def _write_triangles_as_obj_faces(file_like, vertices_at_face):
    """Write triangular patches as OBJ faces.

    OBJ format uses vertex/UVtexture/normal; latter two left blank here.

    Parameters
    ----------
    file_like : file-like
        Opened file-like object to write faces to.
    vertices_at_face : array-like of int, shape *(n_faces, 3)*
        Vertex ids for each triangle.
    """
    for vertex_0, vertex_1, vertex_2 in vertices_at_face:
        print(f"f {vertex_0}// {vertex_1}// {vertex_2}//", file=file_like)


def _write_quads_as_obj_faces(file_, vertices_at_face):
    """Write quad patches as OBJ (triangle) faces.

    Subdivide each patch into two triangles using vertices 0-1-2 (upper left)
    and 2-3-0 (lower right).

    OBJ format uses vertex/UVtexture/normal; latter two left blank here.

    Parameters
    ----------
    file_like : file-like
        Opened file-like object to write faces to.
    vertices_at_face : array-like of int, shape *(n_faces, 4)*
        Vertex ids for each quadrilateral.
    """
    for vertex_0, vertex_1, vertex_2, vertex_3 in vertices_at_face:
        print(f"f {vertex_0}// {vertex_1}// {vertex_2}//", file=file_)
        print(f"f {vertex_2}// {vertex_3}// {vertex_0}//", file=file_)


def _patches_are_triangles(grid):
    """Returns True if patches have 3 nodes, False otherwise."""
    return grid.nodes_at_patch.shape[1] == 3


def _patches_are_quads(grid):
    """Returns True if patches have 4 nodes, False otherwise."""
    return grid.nodes_at_patch.shape[1] == 4


def write_obj(path, grid, field_for_z="topographic__elevation", clobber=False):
    """Write landlab grid and one field to Wavefront OBJ.

    Parameters
    ----------
    path : str, or file-like
        Path to output file.
    grid : Landlab grid object
        Landlab grid object that includes associated values.
    field_for_z : str, optional
        Name of a field to use for the *z*-values of the OBJ file.
    clobber : boolean, optional
        If *path* exists, clobber the existing file, otherwise raise an
        exception.

    Returns
    -------
    str or file-like
        The input path used to write the OBJ file.

    Examples
    --------
    >>> import io
    >>> from landlab import HexModelGrid, RasterModelGrid
    >>> from landlab.io import write_obj

    >>> grid = HexModelGrid((3, 2), spacing=2.0)
    >>> z = grid.add_zeros("topographic__elevation", at="node")
    >>> z[3] = 1.0

    >>> obj_file = write_obj(io.StringIO(), grid)
    >>> print(obj_file.getvalue())
    # landlabgrid
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
    <BLANKLINE>
    """
    if not (_patches_are_triangles(grid) or _patches_are_quads(grid)):
        raise TypeError("grid faces must be triangles or quads")

    z_at_node = grid.at_node[field_for_z]

    if isinstance(path, (str, pathlib.Path)):
        if os.path.exists(path) and not clobber:
            raise ValueError(f"file exists ({path})")

        with open(path, "w") as fp:
            _write_obj_to_filelike(fp, grid, z_at_node)
    else:
        _write_obj_to_filelike(path, grid, z_at_node)

    return path


def _write_obj_to_filelike(file_like, grid, z_at_node):
    file_like.write("# landlabgrid\n")
    file_like.write("#\n")
    file_like.write("g landlabgrid\n")

    _write_obj_vertices(file_like, zip(grid.x_of_node, grid.y_of_node, z_at_node))
    if _patches_are_triangles(grid):
        _write_triangles_as_obj_faces(file_like, grid.nodes_at_patch + 1)
    elif _patches_are_quads(grid):
        _write_quads_as_obj_faces(file_like, grid.nodes_at_patch + 1)

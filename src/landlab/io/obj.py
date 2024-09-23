#! /usr/bin/env python
"""Write (x,y,z) data from a Landlab grid + 1 field to a Wavefront OBJ file.

OBJ functions
+++++++++++++

.. autosummary::

    ~landlab.io.obj.write_obj
"""
import os
import pathlib


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
    >>> from landlab.io.obj import write_obj

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

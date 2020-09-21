#! /usr/bin/env python
"""Write (x,y,z) data from a Landlab grid + 1 field to a Wavefront OBJ file.

OBJ functions
+++++++++++++

.. autosummary::

    ~landlab.io.obj.write_obj
"""

import os


def _write_nodes_as_vertices(file, grid, z):
    """Write node (x,y,z) coordinates as OBJ vertices."""
    for i in range(grid.number_of_nodes):
        file.write(
            "v "
            + str(grid.x_of_node[i])
            + " "
            + str(grid.y_of_node[i])
            + " "
            + str(z[i])
            + "\n"
        )


def _write_triangular_patches_as_faces(file, grid, z):
    """Write triangular patches as OBJ faces.

    OBJ format uses vertex/UVtexture/normal; latter two left blank here."""
    for i in range(len(grid.nodes_at_patch)):
        file.write(
            "f "
            + str(grid.nodes_at_patch[i, 0] + 1)
            + "// "
            + str(grid.nodes_at_patch[i, 1] + 1)
            + "// "
            + str(grid.nodes_at_patch[i, 2] + 1)
            + "//\n"
        )


def _write_quad_patches_as_triangular_faces(file, grid, z):
    """Write quad patches as OBJ (triangle) faces.

    Subdivide each patch into two triangles using vertices 0-1-2 (upper left)
    and 2-3-0 (lower right).

    OBJ format uses vertex/UVtexture/normal; latter two left blank here."""
    for i in range(len(grid.nodes_at_patch)):
        file.write(
            "f "
            + str(grid.nodes_at_patch[i, 0] + 1)
            + "// "
            + str(grid.nodes_at_patch[i, 1] + 1)
            + "// "
            + str(grid.nodes_at_patch[i, 2] + 1)
            + "//\n"
        )
        file.write(
            "f "
            + str(grid.nodes_at_patch[i, 2] + 1)
            + "// "
            + str(grid.nodes_at_patch[i, 3] + 1)
            + "// "
            + str(grid.nodes_at_patch[i, 0] + 1)
            + "//\n"
        )


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
    path : str
        Path to output file.
    grid : Landlab grid object
        Landlab grid object that includes associated values.
    clobber : boolean
        If *path* exists, clobber the existing file, otherwise raise an
        exception.

    Examples
    --------
    >>> import numpy as np
    >>> import os
    >>> import tempfile
    >>> from landlab import HexModelGrid, RasterModelGrid
    >>> from landlab.io.obj import write_obj

    >>> grid = HexModelGrid((3, 2), spacing=2.)
    >>> z = grid.add_zeros("topographic__elevation", at="node")
    >>> z[3] = 1.0
    >>> with tempfile.TemporaryDirectory() as tmpdirname:
    ...     fname = os.path.join(tmpdirname, 'test_tri.obj')
    ...     file = write_obj(fname, grid)
    >>> print(os.path.basename(file))
    test_tri.obj

    >>> grid = RasterModelGrid((3, 3), xy_spacing=2.)
    >>> z = grid.add_zeros("topographic__elevation", at="node")
    >>> z[3] = 1.0
    >>> with tempfile.TemporaryDirectory() as tmpdirname:
    ...     fname = os.path.join(tmpdirname, 'test_quad.obj')
    ...     file = write_obj(fname, grid)
    >>> print(os.path.basename(file))
    test_quad.obj
    """
    if os.path.exists(path) and not clobber:
        raise ValueError("file exists")

    z = grid.at_node[field_for_z]
    with open(path, "w") as f:
        f.write("# landlabgrid\n")
        f.write("#\n")
        f.write("g landlabgrid\n")
        _write_nodes_as_vertices(f, grid, z)
        if _patches_are_triangles(grid):
            _write_triangular_patches_as_faces(f, grid, z)
        elif _patches_are_quads(grid):
            _write_quad_patches_as_triangular_faces(f, grid, z)
        else:
            raise TypeError("grid faces must be triangles or quads")
        f.close()

    return path

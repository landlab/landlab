#!/usr/bin/python
"""
RefinableIcosahedron: class that creates the vertices and triangular
faces of an icosahedron (20-sided polygonal solid) that can be
iteratively refined by replacing each triangle with four smaller
triangles. Designed to provide the starting point for a Landlab
icosphere grid.

Greg Tucker, University of Colorado Boulder, 2023
Adapted from a blog post and code snippet by Andreas Kahler, and
translated into Python (plus function to output in vtk format)
"""

import os
import pathlib

from landlab.io.legacy_vtk import _format_vtk_cells
from landlab.io.legacy_vtk import _format_vtk_header
from landlab.io.legacy_vtk import _format_vtk_points


class RefinableIcosahedron:
    """
    An icosahedron with faces that can be iteratively subdivided.

    Class includes Cartesian coordinates of vertices (12 if not
    subdivided) and IDs of vertices composing each face (20 faces if
    not subdivided). Adapted from a blog post by Andreas Kahler.

    Parameters
    ----------
    radius : float
        Radius for the RefinableIcosahedron, length units (default 1)

    Examples
    --------

    Basic icosahedron

    >>> ico = RefinableIcosahedron()
    >>> len(ico.faces)
    20
    >>> len(ico.vertices)
    12

    Icosphere with two iterations of refinement

    >>> ico.refine_triangles(recursion_level=2)
    >>> len(ico.faces)
    320
    >>> len(ico.vertices)
    162

    Icosahedron with radius != 1

    >>> ico = RefinableIcosahedron(radius=2.0)
    >>> round(ico.vertices[1][0], 2)
    1.05
    >>> round(ico.vertices[0][1], 2)
    1.7
    """

    def __init__(self, radius=1.0):
        """
        Initialize RefinableIcosahedron
        """
        self.radius = radius
        self.vertices = []
        self.faces = []
        self.middle_point_index_cache = {}
        self.vertex_index = -1
        self.create_vertices()
        self.create_faces()

    def add_vertex(self, vtx):
        """
        Add a vertex, scaling its coordinates to fit the given radius,
        and return its index in the vertices list.

        Parameters
        ----------
        vtx : 3-element tuple of float
            x, y, and z coordinates of vertex

        Returns
        -------
        int : index number of the new vertex
        """
        length = (vtx[0] ** 2 + vtx[1] ** 2 + vtx[2] ** 2) ** 0.5
        scale_fac = self.radius / length
        self.vertices.append(
            (vtx[0] * scale_fac, vtx[1] * scale_fac, vtx[2] * scale_fac)
        )
        self.vertex_index += 1
        return self.vertex_index

    def create_vertices(self):
        """
        Create the 12 vertices of an icosahedron.

        Note that the vertex coordinates will be scaled to match
        the given radius.
        """
        t = (1.0 + 5.0**0.5) / 2.0  # this is the famous Golden Mean

        self.add_vertex((-1.0, t, 0.0))
        self.add_vertex((1.0, t, 0.0))
        self.add_vertex((-1.0, -t, 0.0))
        self.add_vertex((1.0, -t, 0.0))

        self.add_vertex((0.0, -1.0, t))
        self.add_vertex((0.0, 1.0, t))
        self.add_vertex((0.0, -1.0, -t))
        self.add_vertex((0.0, 1.0, -t))

        self.add_vertex((t, 0.0, -1.0))
        self.add_vertex((t, 0.0, 1.0))
        self.add_vertex((-t, 0.0, -1.0))
        self.add_vertex((-t, 0.0, 1.0))

    def create_faces(self):
        """
        Create the 20 triangular faces of the icosahedron.

        Faces are stored in a list of 3-element tuples with the vertex IDs.

        Note that "face" here means a triangle on the surface of the object,
        and is different from the Landlab definition of face as the edge
        shared by two neighboring cells.
        """
        # 5 faces around point 0
        self.faces.append((0, 11, 5))
        self.faces.append((0, 5, 1))
        self.faces.append((0, 1, 7))
        self.faces.append((0, 7, 10))
        self.faces.append((0, 10, 11))

        # 5 adjacent faces
        self.faces.append((1, 5, 9))
        self.faces.append((5, 11, 4))
        self.faces.append((11, 10, 2))
        self.faces.append((10, 7, 6))
        self.faces.append((7, 1, 8))

        # 5 faces around point 3
        self.faces.append((3, 9, 4))
        self.faces.append((3, 4, 2))
        self.faces.append((3, 2, 6))
        self.faces.append((3, 6, 8))
        self.faces.append((3, 8, 9))

        # 5 adjacent faces
        self.faces.append((4, 9, 5))
        self.faces.append((2, 4, 11))
        self.faces.append((6, 2, 10))
        self.faces.append((8, 6, 7))
        self.faces.append((9, 8, 1))

    def write_to_vtk(self, path, clobber=False):
        """Save the geometry in a vtk-format text file.

        Note: this function is intended to test the RefinableIcosahedron.
        To write vtk for a Landlab IcosphereGlobalGrid, use
        `landlab.io.legacy_vtk.dump` to capture the full set of geometric
        primitives.

        Parameters
        ----------
        path : str, path-like , or file-like
            Target for output.
        clobber : bool, optional
            Whether to allow overwriting of existing file.

        Returns
        -------
        path : same as input above
            The given output path

        Examples
        --------
        >>> import io
        >>> import os

        >>> ico = RefinableIcosahedron()
        >>> output = ico.write_to_vtk(io.StringIO())
        >>> lines = output.getvalue().splitlines()

        >>> print(lines[0])
        # vtk DataFile Version 2.0

        >>> print(os.linesep.join(lines[5:18]))
        POINTS 12 float
        -0.5257311121191336 0.85065080835204 0.0
        0.5257311121191336 0.85065080835204 0.0
        -0.5257311121191336 -0.85065080835204 0.0
        0.5257311121191336 -0.85065080835204 0.0
        0.0 -0.5257311121191336 0.85065080835204
        0.0 0.5257311121191336 0.85065080835204
        0.0 -0.5257311121191336 -0.85065080835204
        0.0 0.5257311121191336 -0.85065080835204
        0.85065080835204 0.0 -0.5257311121191336
        0.85065080835204 0.0 0.5257311121191336
        -0.85065080835204 0.0 -0.5257311121191336
        -0.85065080835204 0.0 0.5257311121191336

        >>> print(os.linesep.join(lines[18:40]))
        CELLS 20 80
        3 0 11 5
        3 0 5 1
        3 0 1 7
        3 0 7 10
        3 0 10 11
        3 1 5 9
        3 5 11 4
        3 11 10 2
        3 10 7 6
        3 7 1 8
        3 3 9 4
        3 3 4 2
        3 3 2 6
        3 3 6 8
        3 3 8 9
        3 4 9 5
        3 2 4 11
        3 6 2 10
        3 8 6 7
        3 9 8 1

        >>> print(os.linesep.join(lines[41:45]))
        CELL_TYPES 20
        5
        5
        5
        """
        if isinstance(path, (str, pathlib.Path)):
            if os.path.exists(path) and not clobber:
                raise ValueError(f"file exists ({path})")

            with open(path, "w") as fp:
                self._write_vtk_to_filelike(fp)
        else:
            self._write_vtk_to_filelike(path)

        return path

    def _write_vtk_to_filelike(self, file_like):
        """
        Write legacy vtk format to a given file-like object.

        Parameters
        ----------
        file_like : a file-like object (e.g., file pointer, StringIO object)
            The file-like object to write to.

        Returns
        -------
        None
        """
        file_like.write(
            (2 * "\n").join(
                [
                    _format_vtk_header(),
                    _format_vtk_points(self.vertices),
                    _format_vtk_cells(self.faces),
                ]
            )
        )

    def get_middle_point(self, p1, p2):
        """
        Identify and add a new point between two existing points.

        Parameters
        ----------
        p1, p2 : int
            IDs of two points (vertices)

        Returns
        -------
        int : index of the point (vertex) added

        Notes
        -----
        This function is used to refine the icosphere, and is called
        to place new vertices in the middle of the edges of triangles.
        Because each edge is shared by two triangles, we have to avoid
        adding the same point twice. To do this, we use a dictionary
        (middle_point_index_cache) to keep track of points already
        added. Each point has a key made of the two vertex IDs (one
        is bit-shifted, then they are added together). We only add
        a point if it isn't already in the dict.
        """

        # do we already have it?
        key = (min(p1, p2) << 32) + max(p1, p2)
        if key in self.middle_point_index_cache:
            return self.middle_point_index_cache[key]

        # not in cache, so calculate
        point1 = self.vertices[p1]
        point2 = self.vertices[p2]
        middle = (
            (point1[0] + point2[0]) / 2.0,
            (point1[1] + point2[1]) / 2.0,
            (point1[2] + point2[2]) / 2.0,
        )

        i = self.add_vertex(middle)

        # store it and return index
        self.middle_point_index_cache[key] = i
        return i

    def refine_triangles(self, recursion_level=1):
        """
        Subdivide each triangle into four, and add corresponding vertices.

        Parameters
        ----------
        recursion_level : int, optional
            Number of subdivisions to apply (default 1)
        """
        for _ in range(recursion_level):
            faces2 = []
            for tri in self.faces:
                # replace triangle by 4 triangles
                a = self.get_middle_point(tri[0], tri[1])
                b = self.get_middle_point(tri[1], tri[2])
                c = self.get_middle_point(tri[2], tri[0])

                faces2.append((tri[0], a, c))
                faces2.append((tri[1], b, a))
                faces2.append((tri[2], c, b))
                faces2.append((a, b, c))

            self.faces = faces2

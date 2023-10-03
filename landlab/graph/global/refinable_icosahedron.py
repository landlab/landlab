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
        #print("RADIUS", self.radius)
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
        """
        length = (vtx[0] ** 2 + vtx[1] ** 2 + vtx[2] ** 2) ** 0.5
        scale_fac = self.radius / length
        self.vertices.append((vtx[0] * scale_fac, vtx[1] * scale_fac, vtx[2] * scale_fac))
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

        #print("VERTS", self.vertices)

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

    def write_to_vtk(self, filename="icosahedron.vtk"):
        """
        Save the geometry in a vtk-format text file.

        Parameters
        ----------
        filename : str, optional
            Name for output file (defaults to "icosahedron.vtk")
        """
        with open(filename, "w") as f:
            f.write("# vtk DataFile Version 2.0\n")
            f.write("icosahedron\n")
            f.write("ASCII\n")
            f.write("DATASET UNSTRUCTURED_GRID\n")
            f.write("POINTS " + str(len(self.vertices)) + " float\n")
            for vtx in self.vertices:
                f.write(str(vtx[0]) + " " + str(vtx[1]) + " " + str(vtx[2]) + "\n")
            f.write("\n")
            nfaces = len(self.faces)
            f.write("CELLS " + str(nfaces) + " " + str(4 * nfaces) + "\n")
            for face in self.faces:
                f.write(
                    "3 " + str(face[0]) + " " + str(face[1]) + " " + str(face[2]) + "\n"
                )
            f.write("\n")
            f.write("CELL_TYPES " + str(nfaces) + "\n")
            for i in range(nfaces):
                f.write("5\n")
            f.close()

    def dist(self, p1, p2):
        """
        Calculate and return the distance between 3D points p1 and p2.

        Parameters
        ----------
        p1, p2 : tuple, list, or array with 3 float elements
            x, y, z coordinates of each point, m
        """
        dx = p1[0] - p2[0]
        dy = p1[1] - p2[1]
        dz = p1[2] - p2[2]
        return (dx * dx + dy * dy + dz * dz) ** 0.5

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
        is bit shifted, then they are added together). We only add
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
        #print("after refine, verts", self.vertices)
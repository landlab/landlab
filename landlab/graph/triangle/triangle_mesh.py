"""Generates a Delaunay triangulation to be used as a computational mesh.

Given a set of triangulation points (nodes), a bounding box, or the path
to a shapefile, generate a conforming Delaunay triangulation and the
accompanying Voronoi tesselation over the region of interest.

Uses a Python wrapper around Jonathan Shewchuk's Triangle software to
generate the mesh. Documentation for Triangle is located here:

https://www.cs.cmu.edu/~quake/triangle.html

and a list of command line switches that can be passed as opts are here:

https://www.cs.cmu.edu/~quake/triangle.switch.html
"""

import numpy as np
import geopandas as gpd
import shapely
import triangle

from typing import Tuple


class TriangleMesh:
    """Constructs a mesh from a conforming Delaunay triangulation."""

    default_opts = "pqDez"

    def __init__(
        self, points: np.ndarray, segments: list = None, opts: str = default_opts
    ):
        """Initialize the class with a list of triangulation points."""
        self._vertices = np.array(points)
        self._boundary = shapely.LinearRing(self._vertices)

        if segments is None:
            self._segments = self.segment(self._boundary)
        else:
            self._segments = segments

        self.opts = opts

        self.delaunay, self.voronoi = self.triangulate()

    @classmethod
    def from_dims(
        cls,
        shape: Tuple[int, int],
        spacing: Tuple[float, float],
        opts: str = default_opts,
    ):
        """Construct a set of points from shape (nx, ny) and spacing (dx, dy)."""
        nx, ny = shape
        dx, dy = spacing

        n_points = int(nx * ny)
        x_nodes = np.linspace(0, nx * dx, num=nx)
        y_nodes = np.linspace(0, ny * dy, num=ny)

        all_xvals = np.ravel(np.meshgrid(x_nodes, y_nodes)[0])
        all_yvals = np.ravel(np.meshgrid(x_nodes, y_nodes)[1])

        points = np.array([[all_xvals[i], all_yvals[i]] for i in np.arange(n_points)])

        return cls(points, opts=opts)

    @classmethod
    def from_shapefile(cls, path_to_file: str, opts: str = default_opts):
        """Construct a set of points from the boundary of a shapefile."""
        shape = gpd.read_file(path_to_file).geometry
        polygon = shapely.build_area(shape.geometry[0])
        vertices = shapely.get_coordinates(polygon.exterior)

        return cls(vertices, opts=opts)

    def segment(self, curve):
        """Given a LinearRing, return a list of line segments."""
        lines = list(map(shapely.LineString, zip(curve.coords[:-1], curve.coords[1:])))
        segments = []

        for line in lines:
            x1, y1 = line.coords[0]
            x2, y2 = line.coords[1]

            start_vertex = np.argwhere(
                (self._vertices[:, 0] == x1) & (self._vertices[:, 1] == y1)
            )[0]
            end_vertex = np.argwhere(
                (self._vertices[:, 0] == x2) & (self._vertices[:, 1] == y2)
            )[0]

            segments.append([int(start_vertex[0]), int(end_vertex[0])])

        return segments

    def triangulate(self):
        """Call Triangle for this instance's sets of vertices and segments."""
        geometry = {"vertices": self._vertices, "segments": self._segments}

        # We need Triangle to return information about edges
        if "e" not in self.opts:
            self.opts += "e"

        # Python indexes from zero
        if "z" not in self.opts:
            self.opts += "z"

        # Omitting the quality flag will lead to bad meshes
        if "q" not in self.opts:
            raise Warning("Cannot guarantee mesh quality: consider adding 'q' to opts.")

        # Most use cases probably involve Planar Straight Line Graphs
        if "p" not in self.opts:
            raise Warning(
                "If your region is a Planar Straight Line Graph, add 'p' to opts."
            )

        # And, users probably want a conforming Delaunay triangulation
        if "D" not in self.opts:
            raise Warning(
                "If you want a conforming Delaunay triangulation, add 'D' to opts."
            )

        delaunay = triangle.triangulate(geometry, opts=self.opts)

        # We discard ray_origin and ray_direction because by convention the
        # Delaunay vertices form the outermost boundary of the computational grid.
        # Note: this means that number_of_links > number_of_faces.
        points, edges, ray_origin, ray_direction = triangle.voronoi(
            delaunay["vertices"]
        )

        voronoi = {"vertices": points, "edges": edges}

        return delaunay, voronoi

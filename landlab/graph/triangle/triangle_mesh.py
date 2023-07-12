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

from typing import Tuple

class TriangleMesh:
    """Constructs a mesh from a conforming Delaunay triangulation."""

    def __init__(self, points: np.ndarray, segments: list = None):
        """Initialize the class with a list of triangulation points."""
        self._vertices = np.array(points)
        self._boundary = shapely.LinearRing(self._vertices)

        if segments is None:
            self._segments = self.segment(self._boundary)
        else:
            self._segments = segments

    @classmethod
    def from_dims(cls, shape: Tuple[int, int], spacing: Tuple[float, float]):
        """Construct a set of points from shape (nx, ny) and spacing (dx, dy)."""
        nx, ny = shape
        dx, dy = spacing

        n_points = int(nx * ny)
        x_nodes = np.linspace(0, nx * dx, num = nx)
        y_nodes = np.linspace(0, ny * dy, num = ny)

        all_xvals = np.ravel(np.meshgrid(x_nodes, y_nodes)[0])
        all_yvals = np.ravel(np.meshgrid(x_nodes, y_nodes)[1])

        points = np.array([[all_xvals[i], all_yvals[i]] for i in np.arange(n_points)])

        return cls(points)

    @classmethod
    def from_shapefile(cls, path_to_file: str):
        """Construct a set of points from the boundary of a shapefile."""
        shape = gpd.read_file(path_to_file).geometry
        polygon = shapely.build_area(shape.geometry[0])
        vertices = shapely.get_coordinates(polygon.exterior)

        return cls(vertices)

    def segment(self, curve):
        """Given a LinearRing, return a list of line segments."""
        lines = list(map(shapely.LineString, zip(curve.coords[:-1], curve.coords[1:])))
        segments = []

        for line in lines:
            x1, y1 = line.coords[0]
            x2, y2 = line.coords[1]
            
            start_vertex = np.argwhere((self._vertices[:, 0] == x1) & (self._vertices[:, 1] == y1))[0]
            end_vertex = np.argwhere((self._vertices[:, 0] == x2) & (self._vertices[:, 1] == y2))[0]

            segments.append([int(start_vertex[0]), int(end_vertex[0])])

        return segments


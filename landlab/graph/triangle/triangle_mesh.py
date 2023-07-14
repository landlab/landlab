"""Generates a Delaunay triangulation to be used as a computational mesh.

Given a set of triangulation points (nodes), a bounding box, or the path
to a shapefile, generate a conforming Delaunay triangulation and the
accompanying Voronoi tesselation over the region of interest.

Wraps Jonathan Shewchuk's Triangle software to generate the mesh. 
Documentation for Triangle is located here:

https://www.cs.cmu.edu/~quake/triangle.html

and a list of command line switches that can be passed as opts are here:

https://www.cs.cmu.edu/~quake/triangle.switch.html
"""

import numpy as np
import geopandas as gpd
import shapely
from typing import Tuple

class TriangleMesh:
    """Calls Triangle to generate a mesh from a shapefile or array of points."""

    # By default, we want a quality (q) conforming Delaunay triangulation (D)
    # of a polygon (p) with information about edges (e), that indexes from zero (z).
    default_opts = 'pqDez'

    def __init__(self, poly: shapely.Polygon, opts: str = default_opts):
        """Initialize this instance with a Shapely polygon object."""
        self._poly = poly
        self._vertices = shapely.get_coordinates(self._poly)
        self._segments = self._segment(self._poly)
        self._holes = self._identify_holes(self._poly)
        self._opts = opts

        # Dictionaries that are constructed by triangulate()
        self.delaunay = None
        self.voronoi = None

    @classmethod
    def from_shapefile(cls, path_to_file: str, opts: str = default_opts):
        """Initialize this instance with a path to a shapefile."""
        shape = gpd.read_file(path_to_file).geometry
        
        if len(shape) != 1:
            raise TypeError(
                "Shapefile must represent exactly 1 object. \n" +
                "Check that there is only one input geometry and try again."
            )

        if len(shape[0].geoms) != 1:
            raise TypeError(
                "Each shape within the input shapefile must contain exactly 1 geometry."
            )

        polygon = shape[0].geoms[0]

        return cls(polygon, opts)

    @classmethod
    def from_points(
        cls, 
        points: np.ndarray, 
        holes: np.ndarray = None, 
        opts: str = default_opts
    ):
        """Initialize this instance with an array of (x, y) coordinates."""
        polygon = shapely.Polygon(points, holes = holes)
        return cls(polygon, opts = opts)

    def _segment(self, poly: shapely.Polygon) -> np.ndarray:
        """Given a Polygon, construct an array of line segments for exterior and interior rings."""
        boundaries = [poly.exterior]
        interiors = [i for i in poly.interiors]
        segments = []

        boundaries.extend(interiors)

        for curve in boundaries:
            lines = list(map(shapely.LineString, zip(curve.coords[:-1], curve.coords[1:])))

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

        return np.array(segments)

    def _identify_holes(self, shape: shapely.Polygon):
        """Identify interior boundaries within a source Polygon."""
        interiors = [i for i in shape.interiors]
        holes = []

        if len(interiors) > 0:
            for ring in interiors:
                area = shapely.build_area(ring)
                point = area.centroid.xy
                holes.append([point[0][0], point[1][0]])
        else:
            holes = None

        return np.array(holes)
            

    def _write_poly_file(
        self, 
        vertices: np.ndarray, 
        segments: np.ndarray,
        holes: np.ndarray
    ):
        """Write an input .poly file for Triangle."""
        pass

    def _write_mesh_to_file(
        self, 
        vertices: np.ndarray, 
        segments: np.ndarray,
        holes: np.ndarray,
        triangles: np.ndarray
    ):
        """Once a mesh has been created, write the corresponding .node, .ele, and .poly files."""
        pass

    def _read_mesh_files(self, node: str, edge: str, ele: str, v_node: str, v_edge: str):
        """Read output from mesh files."""
        pass

    def triangulate(self, opts: str = default_opts) -> Tuple[dict, dict]:
        """Perform the Delaunay triangulation."""
        delaunay, voronoi = (None, None)

        return delaunay, voronoi

    def refine_mesh(
        self, 
        node_file: str,
        ele_file: str,
        poly_file: str,
        area: np.ndarray
    ) -> Tuple[dict, dict]:
        """Refine the mesh given a new set of maximum area constraints."""
        delaunay, voronoi = (None, None)

        return delaunay, voronoi

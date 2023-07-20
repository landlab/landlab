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

import shutil
import subprocess
import tempfile

import geopandas as gpd
import numpy as np
import pandas as pd
import shapely


class TriangleMesh:
    """Calls Triangle to generate a mesh from a shapefile or array of points."""

    # By default, we want a quality (q) conforming Delaunay triangulation (D)
    # of a polygon (p) with information about edges (e), that indexes from zero (z).
    default_opts = "pqDevz"

    def __init__(self, poly: shapely.Polygon, opts: str = default_opts, timeout=10):
        """Initialize this instance with a Shapely polygon object."""
        self._poly = poly
        self._vertices = shapely.get_coordinates(self._poly)
        self._segments = self._segment(self._poly)
        self._holes = self._identify_holes(self._poly)
        self._opts = opts  # Command-line options to pass to Triangle
        self._timeout = timeout  # How long to let Triangle run before terminating

        # Dictionaries that are constructed by triangulate()
        self.delaunay = None
        self.voronoi = None

    @classmethod
    def from_shapefile(cls, path_to_file: str, opts: str = default_opts, timeout=10):
        """Initialize this instance with a path to a shapefile."""
        shape = gpd.read_file(path_to_file).geometry

        if len(shape) != 1:
            raise TypeError(
                "Shapefile must represent exactly 1 object."
                " Check that there is only one input geometry and try again."
            )

        if len(shape[0].geoms) != 1:
            raise TypeError(
                "Each shape within the input shapefile must contain exactly 1 geometry."
            )

        polygon = shape[0].geoms[0]

        return cls(polygon, opts=opts)

    @classmethod
    def from_points(
        cls,
        points: np.ndarray,
        holes: np.ndarray = None,
        opts: str = default_opts,
        timeout=10,
    ):
        """Initialize this instance with an array of (x, y) coordinates."""
        polygon = shapely.Polygon(points, holes=holes)
        return cls(polygon, opts=opts)

    def _segment(self, poly: shapely.Polygon) -> np.ndarray:
        """Given a Polygon, construct an array of line segments for exterior and
        interior rings.
        """
        boundaries = [poly.exterior]
        interiors = list(poly.interiors)
        segments = []

        boundaries.extend(interiors)

        for curve in boundaries:
            lines = list(
                map(shapely.LineString, zip(curve.coords[:-1], curve.coords[1:]))
            )

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
        interiors = list(shape.interiors)

        holes = []
        for ring in interiors:
            area = shapely.build_area(ring)
            point = area.centroid.xy
            holes.append([point[0][0], point[1][0]])

        return np.array(holes) if len(holes) else None

    def _write_poly_file(
        self, path: str, vertices: np.ndarray, segments: np.ndarray, holes: np.ndarray
    ):
        """Write an input .poly file for Triangle."""
        vertex_header = np.array([vertices.shape[0], 2, 0, 0])[np.newaxis]
        segment_header = np.array([segments.shape[0], 0])[np.newaxis]

        # If there are no holes, the header should just be [0]
        if hasattr(holes, "shape"):
            holes_header = np.array([holes.shape[0]])[np.newaxis]
        else:
            holes_header = np.array([0])

        vertices = np.insert(vertices, 0, np.arange(vertices.shape[0]), axis=1)
        segments = np.insert(segments, 0, np.arange(segments.shape[0]), axis=1)

        # If there are no holes, don't write anything to the .poly file
        if holes_header[0] > 0:
            holes = np.insert(holes, 0, np.arange(holes.shape[0]), axis=1)

        with open(path, "w") as outfile:
            np.savetxt(outfile, vertex_header, fmt="%d")
            np.savetxt(outfile, vertices, fmt="%f")
            np.savetxt(outfile, segment_header, fmt="%d")
            np.savetxt(outfile, segments, fmt="%d")
            np.savetxt(outfile, holes_header, fmt="%d")

            # If there are no holes, there's nothing to write here
            if holes_header[0] > 0:
                np.savetxt(outfile, holes, fmt="%f")

    def _read_mesh_files(
        self, node: str, edge: str, ele: str, v_node: str, v_edge: str
    ) -> tuple[dict, dict]:
        """Read output from mesh files."""
        delaunay = {
            "nodes": pd.read_csv(
                node,
                sep=r"\s+",
                skiprows=1,
                comment="#",
                names=["Node", "x", "y", "BC"],
            ),
            "links": pd.read_csv(
                edge,
                sep=r"\s+",
                skiprows=1,
                comment="#",
                names=["Link", "head", "tail", "BC"],
            ),
            "patches": pd.read_csv(
                ele,
                sep=r"\s+",
                skiprows=1,
                comment="#",
                names=["Patch", "first", "second", "third"],
            ),
        }

        # Triangle writes out rays and edges to the same file,
        # So we need to do some extra work to only keep the Voronoi edges.
        # Read in the data as if everything was defined as a ray,
        faces = pd.read_csv(
            v_edge, sep=r"\s+", skiprows=1, names=["1", "2", "3", "4", "5"], comment="#"
        )[lambda x: x["3"] != -1]
        # then drop any row where the third element ('tail') is undefined.

        # Now we can reshape the array to match the shape we expect from links.
        # Recall that we have discarded any boundary edges from the Voronoi graph.
        faces = faces.drop(["4", "5"], axis=1).rename(
            columns={"1": "Link", "2": "head", "3": "tail"}
        )

        voronoi = {
            "corners": pd.read_csv(
                v_node, sep=r"\s+", skiprows=1, comment="#", names=["Node", "x", "y"]
            ),
            "faces": faces,
        }

        return delaunay, voronoi

    def triangulate(self):
        """Perform the Delaunay triangulation."""

        # -------------------------
        # Check a few items in opts
        # -------------------------
        # We need Triangle to return information about edges
        if "e" not in self._opts:
            self._opts += "e"

        # And information about the Voronoi graph
        if "v" not in self._opts:
            self._opts += "v"

        # Python indexes from zero
        if "z" not in self._opts:
            self._opts += "z"

        # Omitting the quality flag will lead to bad meshes
        if "q" not in self._opts:
            raise Warning("Cannot guarantee mesh quality: consider adding 'q' to opts.")

        # Most use cases probably involve Planar Straight Line Graphs
        if "p" not in self._opts:
            raise Warning(
                "If your region is a Planar Straight Line Graph, add 'p' to opts."
            )

        # And, users probably want a conforming Delaunay triangulation
        if "D" not in self._opts:
            raise Warning(
                "If you want a conforming Delaunay triangulation, add 'D' to opts."
            )

        # --------------------------------
        # Check if Triangle is in the PATH
        # --------------------------------
        path_to_tri = shutil.which("triangle")
        if path_to_tri is None:
            raise OSError(
                "Unable to locate Triangle in PATH. You can install it with:"
                " conda install -c conda-forge triangle"
            )

        # ----------------------------
        # Set up a temporary directory
        # ----------------------------
        with tempfile.TemporaryDirectory() as tmpdir:
            self._write_poly_file(
                tmpdir + "/tri.poly", self._vertices, self._segments, self._holes
            )

            cmd = "triangle"
            input_file = "tri.poly"
            options = "-" + self._opts

            result = subprocess.run(
                [cmd, options, input_file],
                timeout=self._timeout,
                capture_output=True,
                cwd=tmpdir,
            )

            if result.returncode == 0:
                self.delaunay, self.voronoi = self._read_mesh_files(
                    node=tmpdir + "/tri.1.node",
                    edge=tmpdir + "/tri.1.edge",
                    ele=tmpdir + "/tri.1.ele",
                    v_node=tmpdir + "/tri.1.v.node",
                    v_edge=tmpdir + "/tri.1.v.edge",
                )
            else:
                raise OSError(
                    "Triangle failed to generate the mesh, raising the following error:\n"
                    + result.stderr
                )

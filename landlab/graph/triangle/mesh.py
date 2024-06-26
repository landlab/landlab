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

from __future__ import annotations

import pathlib
import shutil
import subprocess
import tempfile

import geopandas as gpd
import numpy as np
import pandas as pd
import shapely


class TriangleMesh:
    """Generate a mesh from a shapefile or array of points.

    Makes a subprocess call to the Triangle command-line interface, which needs
    to be installed separately. The most reliable way to do so is with::

        conda install -c conda-forge triangle

    but users could also install and compile triangle directly.

    Parameters
    ----------
    poly: shapely.Polygon
        The Polygon that serves as the domain for the mesh.
    opts: str
        Command-line arguments to pass to triangle.
    timeout: float
        The amount of time in seconds to let triangle run before terminating.

    See Also
    --------
    TriangleMesh.from_shapefile : Constructor from a shapefile.
    TriangleMesh.from_points : Constructor from an array of (x, y) points.

    Examples
    --------
    >>> poly = shapely.Polygon([[0, 0], [0, 1], [1, 1], [1, 0]])
    >>> tri = TriangleMesh(poly, opts="pqa0.1Devjz", timeout=10)
    >>> tri.triangulate()
    >>> tri.delaunay["nodes"]
        Node     x     y  BC
    0      0  0.00  0.00   1
    1      1  0.00  1.00   1
    2      2  1.00  1.00   1
    3      3  1.00  0.00   1
    4      4  0.50  0.50   0
    5      5  0.00  0.50   1
    6      6  0.50  0.00   1
    7      7  1.00  0.50   1
    8      8  0.50  1.00   1
    9      9  0.75  0.25   0
    10    10  0.25  0.75   0
    11    11  0.75  0.75   0
    12    12  0.25  0.25   0
    >>> tri.voronoi["faces"]
        Link  head  tail
    0      0     0    13
    1      1     0     7
    2      2     0     6
    3      3     1    15
    4      4     1     4
    5      5     1     3
    6      6     2     7
    7      7     2    11
    8      8     2     3
    9      9     3    10
    10    10     4    14
    11    11     4     5
    12    12     5     9
    13    13     5     6
    14    14     6     8
    15    15     7    12
    17    17     8     9
    20    20    10    11
    23    23    12    13
    26    26    14    15
    """

    default_opts = (
        "p"  # Triangulates a Planar Straight Line Graph (.poly file).
        "q"  # Quality mesh generation with no angles smaller than 20 degrees.
        "D"  # Conforming Delaunay.
        "e"  # Outputs (to an .edge file) a list of edges of the triangulation.
        "v"  # Outputs the Voronoi diagram associated with the triangulation.
        "j"  # Jettisons vertices that are not part of the final triangulation.
        "z"  # Numbers all items starting from zero.
    )

    def __init__(self, poly: shapely.Polygon, opts: str = default_opts, timeout=10):
        """Initialize this instance with a Shapely polygon object."""
        self._poly = poly
        self._vertices = shapely.get_coordinates(self._poly)
        self._segments = self._segment(self._poly)
        self._holes = self.identify_holes(self._poly)

        # Command-line options to pass to Triangle
        self._opts = self.validate_options(opts)
        self._timeout = timeout  # How long to let Triangle run before terminating
        self._triangle = self.validate_triangle()

        # Dictionaries that are constructed by triangulate()
        self.delaunay = None
        self.voronoi = None

    @property
    def triangle(self):
        return self._triangle

    @property
    def options(self):
        return self._opts

    @staticmethod
    def validate_options(options):
        for opt in ["e", "v", "z", "j"]:
            options += opt if opt not in options else ""

        # Omitting the quality flag will lead to bad meshes
        if "q" not in options:
            raise Warning("Cannot guarantee mesh quality: consider adding 'q' to opts.")

        # Most use cases probably involve Planar Straight Line Graphs
        if "p" not in options:
            raise Warning(
                "If your region is a Planar Straight Line Graph, add 'p' to opts."
            )

        # And, users probably want a conforming Delaunay triangulation
        if "D" not in options:
            raise Warning(
                "If you want a conforming Delaunay triangulation, add 'D' to opts."
            )

        return options

    @staticmethod
    def validate_triangle(triangle=None):
        triangle = shutil.which("triangle" if triangle is None else triangle)

        if not triangle:
            raise FileNotFoundError(
                "Unable to locate Triangle in PATH. You can install it with:"
                " conda install -c conda-forge triangle"
            )

        subprocess.run([triangle], capture_output=True, check=True)

        return triangle

    @staticmethod
    def read_input_file(path_to_file: str) -> shapely.Polygon:
        """Construct a polygon from an input file."""
        shape = gpd.read_file(path_to_file).geometry

        if len(shape) != 1:
            raise TypeError(
                "Shapefile must represent exactly 1 object."
                " Check that there is only one input geometry and try again."
            )

        if isinstance(shape[0], shapely.MultiPolygon):
            if len(shape[0].geoms) != 1:
                raise TypeError(
                    "Each shape within the input shapefile must contain exactly 1 geometry."
                )

            polygon = shape[0].geoms[0]

        elif isinstance(shape[0], shapely.Polygon):
            polygon = shape[0]

        else:
            raise TypeError(
                "Input geometry is a " + str(type(shape)) + ", but should be a Polygon."
            )

        return polygon

    @staticmethod
    def identify_holes(shape: shapely.Polygon):
        """Identify interior boundaries within a source Polygon."""
        interiors = list(shape.interiors)

        holes = []
        for ring in interiors:
            area = shapely.build_area(ring)
            point = shapely.Point(area.centroid.xy)

            if shape.contains(point):
                invalid = True
                minx, miny, maxx, maxy = ring.bounds

                while invalid:
                    x = np.random.uniform(minx, maxx)
                    y = np.random.uniform(miny, maxy)
                    test_point = shapely.Point(x, y)

                    if not shape.contains(test_point):
                        holes.append([x, y])
                        invalid = False

            else:
                holes.append([point.xy[0][0], point.xy[1][0]])

        return np.array(holes) if len(holes) else None

    @classmethod
    def from_shapefile(cls, path_to_file: str, opts: str = default_opts, timeout=10):
        """Initialize this instance with a path to a shapefile."""
        polygon = cls.read_input_file(path_to_file)
        return cls(polygon, opts=opts, timeout=timeout)

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
        return cls(polygon, opts=opts, timeout=timeout)

    def _segment(self, poly: shapely.Polygon) -> np.ndarray:
        """Given a Polygon, construct an array of line segments for exterior and
        interior rings.
        """
        segments = []

        for ring in poly.interiors:
            lines = list(
                map(shapely.LineString, zip(ring.coords[:-1], ring.coords[1:]))
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

        boundary = list(
            map(
                shapely.LineString,
                zip(poly.exterior.coords[:-1], poly.exterior.coords[1:]),
            )
        )

        for line in boundary:
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

    def _write_poly_file(
        self,
        path: str | pathlib.Path,
        vertices: np.ndarray,
        segments: np.ndarray,
        holes: np.ndarray,
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
            np.savetxt(outfile, vertex_header, fmt="%d", newline="\r\n")
            np.savetxt(outfile, vertices, fmt="%f", newline="\r\n")
            np.savetxt(outfile, segment_header, fmt="%d", newline="\r\n")
            np.savetxt(outfile, segments, fmt="%d", newline="\r\n")
            np.savetxt(outfile, holes_header, fmt="%d", newline="\r\n")

            # If there are no holes, there's nothing to write here
            if holes_header[0] > 0:
                np.savetxt(outfile, holes, fmt="%f")

    def _read_mesh_files(
        self, node: str | pathlib.Path, edge: str, ele: str, v_node: str, v_edge: str
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
        # ----------------------------
        # Set up a temporary directory
        # ----------------------------
        with tempfile.TemporaryDirectory(ignore_cleanup_errors=True) as tmpdir:
            tmp_path = pathlib.Path(tmpdir)

            self._write_poly_file(
                tmp_path / "tri.poly", self._vertices, self._segments, self._holes
            )

            result = subprocess.run(
                [self.triangle, f"-{self.options}", "tri.poly"],
                timeout=self._timeout,
                stdout=subprocess.PIPE,
                stderr=subprocess.STDOUT,
                cwd=tmp_path,
            )

            # if result.returncode == 0:
            try:
                self.delaunay, self.voronoi = self._read_mesh_files(
                    node=tmp_path / "tri.1.node",
                    edge=tmp_path / "tri.1.edge",
                    ele=tmp_path / "tri.1.ele",
                    v_node=tmp_path / "tri.1.v.node",
                    v_edge=tmp_path / "tri.1.v.edge",
                )
            except FileNotFoundError as error:
                raise OSError(
                    "Triangle failed to generate the mesh, raising the following error:\n"
                    + result.stdout.decode()
                ) from error

            # else:
            #     # Triangle sends more informative error messages to stdout
            #     raise OSError(
            #         "Triangle failed to generate the mesh, raising the following error:\n"
            #         + result.stdout.decode()
            #     )

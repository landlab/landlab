import numpy as np
import shapely
import xarray as xr
from shapely.validation import explain_validity

from ..graph import Graph
from ..sort.sort import reverse_one_to_many
from .triangle_mesh import TriangleMesh


class TriangleGraph(Graph):
    """Constructs the Voronoi-Delaunay dual graph.

    We translate Landlab grid elements as follows:

    * nodes: Delaunay vertices
    * links: Delaunay edges
    * patches: Delaunay triangles
    * corners: Voronoi vertices
    * faces: Voronoi edges
    * cells: Voronoi polygons

    Notes
    -----
    By convention, the boundary of the grid will be composed of nodes and links.
    As such, number_of_nodes > number_of_cells and number_of_links > number_of_faces,
    despite the fact that these elements are dual in the Delaunay and Voronoi graphs.

    Parameters
    ----------
    exterior_y_and_x: tuple of array-like
        Coordinates of the exterior nodes, listed with y first, then x.
    holes: array-like
        Coordinates on the boundaries of any interior holes.
    triangle_opts: str
        List of command-line options to pass to Triangle.
    timeout: float
        The amount of time in seconds to let Triangle run before terminating.

    Examples
    --------
    >>> ys = [0, 0, 10, 10]
    >>> xs = [0, 10, 10, 0]
    >>> graph = TriangleGraph(np.array([ys, xs]), triangle_opts="pqa1Djevz")
    >>> graph.number_of_nodes
    89
    >>> graph.number_of_links
    232
    >>> graph.nodes_at_link[:3]
    array([[41, 43],
           [43, 15],
           [15, 41]])
    """

    def __init__(
        self,
        exterior_y_and_x: tuple[np.ndarray, np.ndarray],
        holes: np.ndarray = None,
        triangle_opts: str = "",
        timeout: float = 10,
    ):
        """Initialize this instance with dicts of Delaunay and Voronoi geometries."""
        ys, xs = exterior_y_and_x
        exterior_points = np.column_stack([xs, ys])
        polygon = shapely.Polygon(exterior_points, holes=holes)

        if polygon.is_valid is False:
            raise ValueError(
                "Shapely considers the input geometry invalid for the following reasons:\n"
                + str(explain_validity(polygon))
            )

        mesh_generator = TriangleMesh(
            polygon,
            opts=triangle_opts,
            timeout=timeout,
        )

        mesh_generator.triangulate()

        self._delaunay = mesh_generator.delaunay
        self._voronoi = mesh_generator.voronoi

        for required in ["nodes", "links", "patches"]:
            if required not in self._delaunay:
                raise ValueError(f"Missing {required!r} in Delaunay dictionary.")

        for required in ["corners", "faces"]:
            if required not in self._voronoi:
                raise ValueError(f"Missing {required!r} in Voronoi dictionary.")

        self._mesh = xr.Dataset(
            {
                "node": xr.DataArray(
                    data=np.arange(self._delaunay["nodes"].shape[0]),
                    coords={
                        "x_of_node": xr.DataArray(
                            self._delaunay["nodes"]["x"], dims=("node",)
                        ),
                        "y_of_node": xr.DataArray(
                            self._delaunay["nodes"]["y"], dims=("node",)
                        ),
                    },
                    dims=("node",),
                ),
                "corner": xr.DataArray(
                    data=np.arange(self._voronoi["corners"].shape[0]),
                    coords={
                        "x_of_corner": xr.DataArray(
                            self._voronoi["corners"]["x"], dims=("corner",)
                        ),
                        "y_of_corner": xr.DataArray(
                            self._voronoi["corners"]["y"], dims=("corner",)
                        ),
                    },
                    dims=("corner",),
                ),
            }
        )

        self._mesh.update(
            {
                "nodes_at_link": xr.DataArray(
                    np.dstack(
                        [
                            self._delaunay["links"]["head"].values,
                            self._delaunay["links"]["tail"].values,
                        ]
                    )[0],
                    dims=("link", "Two"),
                )
            }
        )

        nodes_at_cell, cells_at_node = self._number_cells()
        corners_at_cell = self._get_corners_at_cell(nodes_at_cell, cells_at_node)
        nodes_at_patch, links_at_patch = self._get_nodes_and_links_at_patch()
        nodes_at_face = self._get_nodes_at_face()
        faces_at_cell = self._get_faces_at_cell(corners_at_cell)

        self._mesh.update(
            {
                "nodes_at_patch": xr.DataArray(nodes_at_patch, dims=("patch", "Three")),
                "corners_at_face": xr.DataArray(
                    np.dstack(
                        [
                            self._voronoi["faces"]["head"].values,
                            self._voronoi["faces"]["tail"].values,
                        ]
                    )[0],
                    dims=("face", "Two"),
                ),
                "corners_at_cell": xr.DataArray(
                    corners_at_cell, dims=("cell", "max_corners_per_cell")
                ),
                "n_corners_at_cell": xr.DataArray(
                    [len(corners) for corners in corners_at_cell], dims=("cell",)
                ),
                "nodes_at_face": xr.DataArray(nodes_at_face, dims=("face", "Two")),
                "cell_at_node": xr.DataArray(cells_at_node, dims=("node",)),
                "links_at_patch": xr.DataArray(links_at_patch, dims=("patch", "Three")),
                "node_at_cell": xr.DataArray(nodes_at_cell, dims=("cell",)),
                "faces_at_cell": xr.DataArray(
                    faces_at_cell, dims=("cell", "max_faces_per_cell")
                ),
            }
        )

        self.perimeter_nodes = np.where(np.array(self._delaunay["nodes"])[:, 3] == 1)[0]

        Graph.__init__(
            self,
            (self._mesh.y_of_node, self._mesh.x_of_node),
            links=self._mesh.nodes_at_link,
            patches=self._mesh.links_at_patch,
            sort=False,
        )
        dual_graph = Graph.__init__(
            self,
            (self._mesh.y_of_corner, self._mesh.x_of_corner),
            links=self._mesh.corners_at_face,
            patches=self._mesh.faces_at_cell,
            sort=False,
        )

        self.merge(
            dual_graph,
            node_at_cell=self._mesh.node_at_cell,
            nodes_at_face=self._mesh.nodes_at_face,
        )

    @classmethod
    def from_shapefile(
        cls, path_to_file: str, triangle_opts: str = "", timeout: float = 10
    ):
        """Initialize a TriangleGraph from an input file."""
        polygon = TriangleMesh.read_input_file(path_to_file)
        nodes_y = np.array(polygon.exterior.xy[1])
        nodes_x = np.array(polygon.exterior.xy[0])
        holes = polygon.interiors

        return cls(
            (nodes_y, nodes_x),
            holes=holes,
            triangle_opts=triangle_opts,
            timeout=timeout,
        )

    def _number_cells(self) -> tuple[np.ndarray, np.ndarray]:
        """Map between grid cells and nodes."""
        nodes_at_cell = self._delaunay["nodes"]["Node"][
            self._delaunay["nodes"]["BC"] == 0
        ].values
        cells_at_node = np.array(
            [
                np.argwhere(nodes_at_cell == i)[0][0] if i in nodes_at_cell else -1
                for i in np.arange(self.number_of_nodes)
            ]
        )

        return nodes_at_cell, cells_at_node

    def _get_corners_at_cell(self, nodes_at_cell, cells_at_node) -> np.ndarray:
        """Construct an array of size (n_cells, max_corners_per_cell) from the Voronoi graph."""
        patches = np.array(self._delaunay["patches"])

        corners = [[] for _ in range(len(nodes_at_cell))]
        max_corners_per_cell = 0

        for node, cell in enumerate(cells_at_node):
            if cell != -1:
                triangles = (patches[:, 1:] == node).any(axis=1).nonzero()[0]
                corners[cell].extend(triangles)

                if len(triangles) > max_corners_per_cell:
                    max_corners_per_cell = len(triangles)

        corners_at_cell = np.array(
            [
                np.pad(
                    cell_corners,
                    (0, max_corners_per_cell - len(cell_corners)),
                    "constant",
                    constant_values=-1,
                )
                for cell_corners in corners
            ],
            dtype=int,
        )

        if max_corners_per_cell == 0:
            raise ValueError(
                "Triangle failed to generate cells.\n"
                "This is most likely to occur when all nodes are on boundaries,\n"
                "and may be fixed by making a finer mesh,\n"
                "e.g., by setting a smaller value of 'a' in opts."
            )

        return corners_at_cell

    def _get_nodes_at_face(self) -> np.ndarray:
        """Construct an array of size (n_faces, 2) from the Voronoi graph."""
        nodes_at_link = np.array(self._delaunay["links"])
        nodes_at_face = np.full((self._voronoi["faces"].shape[0], 2), -1)
        nodes_at_face[:] = nodes_at_link[: nodes_at_face.shape[0], 1:3]

        return nodes_at_face

    def _get_faces_at_cell(self, corners_at_cell) -> np.ndarray:
        """Construct an array of size (n_cells, max_faces_per_cell) from the Delaunay graph."""
        n_cells, max_corners_per_cell = corners_at_cell.shape
        faces_at_corner = reverse_one_to_many(np.array(self._voronoi["faces"])[:, 1:])

        faces = [[] for _ in range(n_cells)]
        for cell in range(n_cells):
            possible_faces = []

            for corner in corners_at_cell[cell]:
                if corner != -1:
                    possible_faces += list(
                        faces_at_corner[corner][faces_at_corner[corner] != -1]
                    )

            unique, count = np.unique(possible_faces, return_counts=True)
            faces[cell].extend(unique[count > 1])

        faces_at_cell = np.array(
            [
                np.pad(
                    cell_faces,
                    (0, max_corners_per_cell - len(cell_faces)),
                    "constant",
                    constant_values=-1,
                )
                for cell_faces in faces
            ],
            dtype=int,
        )

        return faces_at_cell

    def _get_nodes_and_links_at_patch(self) -> np.ndarray:
        """From the Delaunay graph, identify the nodes and links adjacent to each patch."""
        links_at_patch = np.full((self._delaunay["patches"].shape[0], 3), -1)
        links_at_node = reverse_one_to_many(np.array(self._delaunay["links"])[:, 1:3])
        nodes_at_patch = np.dstack(
            [
                self._delaunay["patches"]["first"],
                self._delaunay["patches"]["second"],
                self._delaunay["patches"]["third"],
            ]
        )[0]

        for patch in np.arange(links_at_patch.shape[0]):
            possible_links = []

            for node in nodes_at_patch[patch]:
                possible_links += list(links_at_node[node][links_at_node[node] != -1])

            unique, count = np.unique(possible_links, return_counts=True)
            links_at_patch[patch, :] = unique[count > 1]

        return nodes_at_patch, links_at_patch

    @property
    def number_of_nodes(self):
        return self._mesh.dims["node"]

    @property
    def number_of_links(self):
        return self._mesh.dims["link"]

    @property
    def number_of_patches(self):
        return self._mesh.dims["patch"]

    @property
    def number_of_corners(self):
        return self._mesh.dims["corner"]

    @property
    def number_of_faces(self):
        return self._mesh.dims["face"]

    @property
    def number_of_cells(self):
        return self._mesh.dims["cell"]

    @property
    def x_of_node(self):
        return self._mesh["x_of_node"].values

    @property
    def y_of_node(self):
        return self._mesh["y_of_node"].values

    @property
    def x_of_corner(self):
        return self._mesh["x_of_corner"].values

    @property
    def y_of_corner(self):
        return self._mesh["y_of_corner"].values

    @property
    def nodes_at_patch(self):
        return self._mesh["nodes_at_patch"].values

    @property
    def nodes_at_link(self):
        return self._mesh["nodes_at_link"].values

    @property
    def nodes_at_face(self):
        return self._mesh["nodes_at_face"].values

    @property
    def corners_at_face(self):
        return self._mesh["corners_at_face"].values

    @property
    def corners_at_cell(self):
        return self._mesh["corners_at_cell"].values

    @property
    def n_corners_at_cell(self):
        return self._mesh["n_corners_at_cell"].values

    @property
    def cell_at_node(self):
        return self._mesh["cell_at_node"].values

    @property
    def links_at_patch(self):
        return self._mesh["links_at_patch"].values

    @property
    def node_at_cell(self):
        return self._mesh["node_at_cell"].values

    @property
    def faces_at_cell(self):
        return self._mesh["faces_at_cell"].values

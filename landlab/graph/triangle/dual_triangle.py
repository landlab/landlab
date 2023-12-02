import numpy as np
import shapely
from shapely.validation import explain_validity

from landlab.graph.dual import DualGraph
from landlab.graph.graph import Graph
from landlab.graph.sort.sort import reverse_one_to_many, sort_links_at_patch
from landlab.graph.triangle.triangle_mesh import TriangleMesh


class TriGraph(Graph):
    def __init__(
        self,
        exterior_y_and_x: tuple[np.ndarray, np.ndarray],
        holes: np.ndarray = None,
        triangle_opts: str = "",
        timeout: float = 10,
        sort: bool = False,
    ):
        polygon = shapely.Polygon(
            zip(exterior_y_and_x[1], exterior_y_and_x[0]), holes=holes
        )

        if not polygon.is_valid:
            raise ValueError(
                "Shapely considers the input geometry invalid for the following"
                f" reasons:\n{explain_validity(polygon)}"
            )

        mesh_generator = TriangleMesh(polygon, opts=triangle_opts, timeout=timeout)

        mesh_generator.triangulate()
        self._delaunay = mesh_generator.delaunay

        nodes_at_link = np.ascontiguousarray(
            self._delaunay["links"][["head", "tail"]].values
        )
        _, links_at_patch = self._get_nodes_and_links_at_patch()

        self.perimeter_nodes = np.flatnonzero(self._delaunay["nodes"]["BC"] == 1)

        Graph.__init__(
            self,
            (self._delaunay["nodes"]["x"], self._delaunay["nodes"]["y"]),
            links=nodes_at_link,
            patches=links_at_patch,
            sort=sort,
        )

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


class DualTriGraph(DualGraph, TriGraph):
    def __init__(
        self,
        exterior_y_and_x: tuple[np.ndarray, np.ndarray],
        holes: np.ndarray = None,
        triangle_opts: str = "",
        timeout: float = 10,
        sort: bool = False,
    ):
        polygon = shapely.Polygon(
            zip(exterior_y_and_x[1], exterior_y_and_x[0]), holes=holes
        )

        if not polygon.is_valid:
            raise ValueError(
                "Shapely considers the input geometry invalid for the following reasons:\n"
                + str(explain_validity(polygon))
            )

        mesh_generator = TriangleMesh(polygon, opts=triangle_opts, timeout=timeout)

        mesh_generator.triangulate()
        self._delaunay = mesh_generator.delaunay

        nodes_at_link = np.ascontiguousarray(
            self._delaunay["links"][["head", "tail"]].values
        )
        _, links_at_patch = self._get_nodes_and_links_at_patch()

        self.perimeter_nodes = np.flatnonzero(self._delaunay["nodes"]["BC"] == 1)

        Graph.__init__(
            self,
            (self._delaunay["nodes"]["y"], self._delaunay["nodes"]["x"]),
            links=nodes_at_link,
            patches=links_at_patch,
            sort=False,
        )

        self._voronoi = mesh_generator.voronoi
        corners_at_face = self._voronoi["faces"][["head", "tail"]].values

        node_at_cell, cells_at_node = self._number_cells()
        nodes_at_face = self._get_nodes_at_face()
        corners_at_cell = self._get_corners_at_cell(node_at_cell, cells_at_node)
        faces_at_cell = self._get_faces_at_cell(corners_at_cell)

        xy_of_corner = np.c_[
            (
                self._voronoi["corners"]["y"],
                self._voronoi["corners"]["x"],
            )
        ]
        sort_links_at_patch(faces_at_cell, corners_at_face, xy_of_corner)

        dual_graph = Graph(
            (
                np.ascontiguousarray(self._voronoi["corners"]["y"]),
                np.ascontiguousarray(self._voronoi["corners"]["x"]),
            ),
            links=np.ascontiguousarray(corners_at_face),
            patches=np.ascontiguousarray(faces_at_cell),
            sort=False,
        )

        self.merge(
            dual_graph,
            node_at_cell=np.ascontiguousarray(node_at_cell),
            nodes_at_face=np.ascontiguousarray(nodes_at_face),
        )

        if sort:
            self.sort()

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

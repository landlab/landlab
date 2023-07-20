import numpy as np
import xarray as xr


class TriangleGraph:
    """Constructs the Voronoi-Delaunay dual graph.

    We translate Landlab grid elements as follows:
    * nodes: Delaunay vertices
    * links: Delaunay edges
    * patches: Delaunay triangles
    * corners: Voronoi vertices
    * faces: Voronoi edges
    * cells: Voronoi polygons

    Note: by convention, the boundary of the grid will be composed of nodes and links.
        As such, number_of_nodes > number_of_cells and number_of_links > number_of_faces,
        despite the fact that these elements are dual in the Delaunay and Voronoi graphs.
    """

    def __init__(self, delaunay: dict, voronoi: dict):
        """Initialize this instance with dicts of Delaunay and Voronoi geometries."""

        for required in ["nodes", "links", "patches"]:
            if required not in delaunay:
                raise ValueError(f"Missing {required!r} in Delaunay dictionary.")

        for required in ["corners", "faces"]:
            if required not in voronoi:
                raise ValueError(f"Missing {required!r} in Voronoi dictionary.")

        self._delaunay = delaunay
        self._voronoi = voronoi

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

        corners_at_cell = self._get_corners_at_cell()

        self._mesh.update(
            {
                "nodes_at_patch": xr.DataArray(
                    self._get_nodes_and_links_at_patch()[0], dims=("patch", "Three")
                ),
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
                "nodes_at_face": xr.DataArray(
                    self._get_nodes_at_face(), dims=("face", "Two")
                ),
                "cell_at_node": xr.DataArray(self._number_cells()[1], dims=("node",)),
                "links_at_patch": xr.DataArray(
                    self._get_nodes_and_links_at_patch()[1], dims=("patch", "Three")
                ),
                "node_at_cell": xr.DataArray(self._number_cells()[0], dims=("cell",)),
                "faces_at_cell": xr.DataArray(
                    self._get_faces_at_cell(), dims=("cell", "max_faces_per_cell")
                ),
            }
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

    def _get_corners_at_cell(self) -> np.ndarray:
        """Construct an array of size (n_cells, max_corners_at_cell) from the Voronoi graph."""
        nodes_at_cell, cells_at_node = self._number_cells()
        max_corners_per_cell = 0

        corners = [[] for _ in nodes_at_cell]
        for node, cell in enumerate(cells_at_node):
            if cell != -1:
                triangles = np.where(
                    np.isin(
                        self._delaunay["patches"].loc[
                            :, self._delaunay["patches"].columns != "Patch"
                        ],
                        node,
                    )
                )[0]
                corners[cell] = triangles

                if len(triangles) > max_corners_per_cell:
                    max_corners_per_cell = len(triangles)

        corners_at_cell = np.full((len(nodes_at_cell), max_corners_per_cell), -1)
        for cell in range(nodes_at_cell.shape[0]):
            corners_at_cell[cell, : len(corners[cell])] = corners[cell]

        return corners_at_cell

    def _get_nodes_at_face(self) -> np.ndarray:
        """Construct an array of size (n_faces, 2) from the Voronoi graph."""
        nodes_at_face = np.full((len(self._voronoi["faces"]), 2), -1)

        for face in np.arange(self._voronoi["faces"].shape[0]):
            if face < self._delaunay["links"].shape[0]:
                nodes_at_face[face] = (
                    self._delaunay["links"].iloc[face]["head"],
                    self._delaunay["links"].iloc[face]["tail"],
                )

        return nodes_at_face

    def _get_faces_at_cell(self) -> np.ndarray:
        """Construct an array of size (n_cells, max_faces_per_cell) from the Delaunay graph."""
        corners_at_cell = self._get_corners_at_cell()

        faces = [[] for _ in range(corners_at_cell.shape[0])]
        max_faces_at_cell = 0

        for cell, corners in enumerate(corners_at_cell):
            possible_faces = []

            for corner in corners:
                possible_faces += list(
                    np.where(np.isin(self._voronoi["faces"]["head"], corner))[0]
                )
                possible_faces += list(
                    np.where(np.isin(self._voronoi["faces"]["tail"], corner))[0]
                )

            unique, count = np.unique(possible_faces, return_counts=True)
            faces[cell] = unique[count > 1]

            if len(unique[count > 1]) > max_faces_at_cell:
                max_faces_at_cell = len(unique[count > 1])

        faces_at_cell = np.full((corners_at_cell.shape[0], max_faces_at_cell), -1)
        for cell in range(faces_at_cell.shape[0]):
            faces_at_cell[cell, : len(faces[cell])] = faces[cell]

        return faces_at_cell

    def _get_nodes_and_links_at_patch(self) -> np.ndarray:
        """From the Delaunay graph, identify the nodes and links adjacent to each patch."""
        links_at_patch = np.full((self._delaunay["patches"].shape[0], 3), -1)
        nodes_at_patch = np.dstack(
            [
                self._delaunay["patches"]["first"],
                self._delaunay["patches"]["second"],
                self._delaunay["patches"]["third"],
            ]
        )[0]

        for patch in np.arange(links_at_patch.shape[0]):
            adjacent_nodes = nodes_at_patch[patch]
            possible_links = []

            for node in adjacent_nodes:
                possible_links += list(
                    np.where(np.isin(self._delaunay["links"]["head"], node))[0]
                )
                possible_links += list(
                    np.where(np.isin(self._delaunay["links"]["tail"], node))[0]
                )

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

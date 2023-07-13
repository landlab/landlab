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

        for required in [
            'vertices', 
            'vertex_markers', 
            'edges',
            'edge_markers',
            'triangles'
        ]:
            if required not in delaunay.keys():
                raise ValueError("Missing '" + str(required) + "' in Delaunay dictionary.")

        for required in [
            'vertices',
            'edges'
        ]:
            if required not in voronoi.keys():
                raise ValueError("Missing '" + str(required) + "' in Voronoi dictionary.")

        self._delaunay = delaunay
        self._voronoi = voronoi

        self._mesh = xr.Dataset(
            {
                "node": xr.DataArray(
                    data = np.arange(len(self._delaunay['vertices'])),
                    coords = {
                        "x_of_node": xr.DataArray(self._delaunay['vertices'][:, 0], dims = ("node",)),
                        "y_of_node": xr.DataArray(self._delaunay['vertices'][:, 1], dims = ("node",)),
                    },
                    dims = ("node",),
                ),
                "corner": xr.DataArray(
                    data = np.arange(len(self._voronoi['vertices'])),
                    coords = {
                        "x_of_corner": xr.DataArray(self._voronoi['vertices'][:, 0], dims = ("corner",)),
                        "y_of_corner": xr.DataArray(self._voronoi['vertices'][:, 1], dims = ("corner",)),
                    },
                    dims = ("corner",),
                )
            }
        )

        self._mesh.update(
            {
                "nodes_at_link": xr.DataArray(
                    self._delaunay['edges'], dims = ("link", "Two")
                ),
                "nodes_at_patch": xr.DataArray(
                    self._delaunay['triangles'], dims = ("patch", "Three")
                ),
                "corners_at_face": xr.DataArray(
                    self._voronoi['edges'], dims = ("face", "Two")
                ),
                # "corners_at_cell": xr.DataArray(
                #     corners_at_cell, dims = ("cell", "max_corners_per_cell")
                # ), 
                # "n_corners_at_cell": xr.DataArray(
                #     [len(cell) for cell in corners_at_cell], dims = ("cell",)
                # ),
                "nodes_at_face": xr.DataArray(
                    self._get_nodes_at_face(), dims = ("face", "Two")
                ),
                # "cell_at_node": xr.DataArray(
                #     np.arange(self.number_of_nodes), dims = ("node",)
                # ),
                "links_at_patch": xr.DataArray(
                    self._delaunay['triangles'], dims = ("patch", "Three")
                ),
                # "node_at_cell": xr.DataArray(
                #     np.arange(self.number_of_nodes), dims = ("cell",)
                # ),
                # "faces_at_cell": xr.DataArray(
                #     self._get_faces_at_cell(), dims = ("cell", "max_faces_per_cell")
                # )
            }
        )

    def _get_corners_at_cell(self) -> np.ndarray:
        """Construct an array of size (n_cells, max_corners_at_cell) from the Voronoi graph."""
        cells = [
            idx for idx in range(len(self._delaunay['vertex_markers']))
            if self._delaunay['vertex_markers'][idx] == 0
        ]
        max_corners_per_cell = 0

        corners = {cell: [] for cell in cells}

        for cell in cells:
            triangles = np.where(np.isin(self._delaunay['triangles'], cell))[0]
            corners[cell] = triangles

            if len(triangles) > max_corners_per_cell:
                max_corners_per_cell = len(triangles)

        # TODO Do we want to keep cell j as the dual of node j, or reindex cells from [0, N] ?
        
        return None

    def _get_nodes_at_face(self) -> np.ndarray:
        """Construct an array of size (n_faces, 2) from the Voronoi graph."""
        nodes_at_face = np.full((len(self._voronoi['edges']), 2), -1)

        for face in np.arange(len(self._voronoi['edges'])):
            if face < len(self._delaunay['edges']):
                nodes_at_face[face] = self._delaunay['edges'][face]
            else:
                pass

        return nodes_at_face

    def _get_faces_at_cell(self) -> np.ndarray:
        """Construct an array of size (n_cells, max_faces_per_cell) from the Delaunay graph."""
        
        # TODO see comment in _get_corners_at_cell
        pass        
            
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
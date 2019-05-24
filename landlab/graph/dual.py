"""Define a graph of nodes-links-patches and its dual.

This class should not be used directly. Instead, it should be used as a base
class when defining other types of graphs.
"""
import numpy as np
import xarray as xr

from ..core.utils import as_id_array
from ..utils.decorators import read_only_array, store_result_in_grid
from .graph import Graph, find_perimeter_nodes
from .sort.sort import reverse_one_to_one


def _sort_dual_graph(graph):
    from .sort.sort import reindex_by_xy
    from .sort.ext.remap_element import remap_graph_element

    sorted_dual = reindex_by_xy(graph._dual)
    sorted = reindex_by_xy(graph)

    node_at_cell = graph.ds["node_at_cell"].values
    node_at_cell[:] = node_at_cell[sorted_dual[2]]
    remap_graph_element(graph.node_at_cell, as_id_array(np.argsort(sorted[0])))


def update_node_at_cell(ugrid, node_at_cell):
    node_at_cell = xr.DataArray(
        data=as_id_array(node_at_cell),
        dims=("cell",),
        attrs={
            "cf_role": "cell_node_connectivity",
            "long_name": "nodes centered at cells",
            "start_index": 0,
        },
    )
    ugrid.update({"node_at_cell": node_at_cell})


def update_nodes_at_face(ugrid, nodes_at_face):
    nodes_at_face = xr.DataArray(
        data=as_id_array(nodes_at_face),
        dims=("face", "Two"),
        attrs={
            "cf_role": "face_node_connectivity",
            "long_name": "nodes on either side of a face",
            "start_index": 0,
        },
    )
    ugrid.update({"nodes_at_face": nodes_at_face})


class DualGraph(Graph):
    def __init__(self, **kwds):
        node_at_cell = kwds.pop("node_at_cell", None)
        nodes_at_face = kwds.pop("nodes_at_face", None)

        update_node_at_cell(self.ds, node_at_cell)
        update_nodes_at_face(self.ds, nodes_at_face)

        rename = {
            "mesh": "dual",
            "node": "corner",
            "link": "face",
            "patch": "cell",
            "x_of_node": "x_of_corner",
            "y_of_node": "y_of_corner",
            "nodes_at_link": "corners_at_face",
            "links_at_patch": "faces_at_cell",
            "max_patch_links": "max_cell_faces",
        }
        self._ds = xr.merge([self._ds, self._dual.ds.rename(rename)])

        self._origin = (0.0, 0.0)

        self._frozen = False
        self.freeze()

        if kwds.get("sort", True):
            self.sort()

    def sort(self):
        from .sort.ext.remap_element import remap_graph_element

        sorted_nodes, sorted_links, sorted_patches = Graph.sort(self)
        sorted_corners, sorted_faces, sorted_cells = self.dual.sort()

        with self.thawed():
            self.node_at_cell[:] = self.node_at_cell[sorted_cells]
            self.nodes_at_face[:] = self.nodes_at_face[sorted_faces]

            remap_graph_element(
                as_id_array(self.node_at_cell), as_id_array(np.argsort(sorted_nodes))
            )
            remap_graph_element(
                as_id_array(self.nodes_at_face).reshape((-1,)),
                as_id_array(np.argsort(sorted_nodes)),
            )

    def freeze(self):
        Graph.freeze(self)
        if hasattr(self, "dual"):
            self.dual.freeze()

    def thaw(self):
        Graph.thaw(self)
        if hasattr(self, "dual"):
            self.dual.thaw()

    @property
    def dual(self):
        return self._dual

    @property
    def node_at_cell(self):
        return self.ds["node_at_cell"].values

    @property
    def nodes_at_face(self):
        return self.ds["nodes_at_face"].values

    @property
    def cell_at_node(self):
        try:
            return self._cell_at_node
        except AttributeError:
            self._cell_at_node = reverse_one_to_one(
                self.node_at_cell, minlength=self.number_of_nodes
            )
            return self._cell_at_node

    @property
    def link_at_face(self):
        try:
            return self._link_at_face
        except AttributeError:
            return self._create_link_at_face()

    def _create_link_at_face(self):

        link_at_nodes = {}
        for link, pair in enumerate(self.nodes_at_link):
            # pair.sort()
            link_at_nodes[tuple(np.sort(pair))] = link

        link_at_face = np.full((self.number_of_faces,), -1, dtype=int)
        # for face, pair in enumerate(self._nodes_at_face):
        for face, pair in enumerate(self.nodes_at_face):
            # pair.sort()
            link_at_face[face] = link_at_nodes[tuple(np.sort(pair))]
        self._link_at_face = link_at_face

        return self._link_at_face

    @property
    def face_at_link(self):
        try:
            return self._face_at_link
        except AttributeError:
            self._face_at_link = reverse_one_to_one(
                self.link_at_face, minlength=self.number_of_links
            )
            return self._face_at_link

    @property
    def x_of_corner(self):
        return self._dual.x_of_node

    @property
    def y_of_corner(self):
        return self._dual.y_of_node

    @property
    def xy_of_corner(self):
        return self._dual.xy_of_node

    @property
    def xy_of_face(self):
        return self._dual.xy_of_link

    @property
    def xy_of_cell(self):
        return self._dual.xy_of_patch

    @property
    def corners(self):
        return self._dual.nodes

    @property
    def number_of_corners(self):
        return self._dual.number_of_nodes

    @property
    @store_result_in_grid()
    @read_only_array
    def perimeter_corners(self):
        return find_perimeter_nodes(self.dual)

    @property
    def corners_at_face(self):
        return self._dual.nodes_at_link

    @property
    def corner_at_face_tail(self):
        return self._dual.node_at_link_tail

    @property
    def corner_at_face_head(self):
        return self._dual.node_at_link_head

    @property
    def number_of_faces(self):
        return self._dual.number_of_links

    @property
    def faces_at_cell(self):
        return self._dual.links_at_patch

    @property
    def corners_at_cell(self):
        return self._dual.nodes_at_patch

    @property
    def number_of_cells(self):
        return self._dual.number_of_patches

    @property
    def faces_at_corner(self):
        return self._dual.links_at_node

    @property
    def face_dirs_at_corner(self):
        return self._dual.link_dirs_at_node

    @property
    def cells_at_face(self):
        return self._dual.patches_at_link

    @property
    def cells_at_corner(self):
        return self._dual.patches_at_node

    @property
    def width_of_face(self):
        return self._dual.length_of_link

    @property
    def length_of_face(self):
        return self._dual.length_of_link

    @property
    def area_of_cell(self):
        return self._dual.area_of_patch

    @property
    def adjacent_corners_at_corner(self):
        return self._dual.adjacent_nodes_at_node

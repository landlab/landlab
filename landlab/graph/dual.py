"""Define a graph of nodes-links-patches and its dual.

This class should not be used directly. Instead, it should be used as a base
class when defining other types of graphs.
"""
import numpy as np

from ..core.utils import as_id_array
from ..utils.decorators import store_result_in_grid, read_only_array
from .graph import Graph, find_perimeter_nodes
from .sort.sort import reverse_one_to_one


def _sort_dual_graph(graph):
    from .sort.sort import reindex_by_xy
    from .sort.ext.remap_element import remap_graph_element

    sorted_dual = reindex_by_xy(graph._dual)
    sorted = reindex_by_xy(graph)

    graph._node_at_cell = graph._node_at_cell[sorted_dual[2]]
    remap_graph_element(graph._node_at_cell,
                        as_id_array(np.argsort(sorted[0])))


class DualGraph(object):
    def __init__(self, *args, **kwds):
        node_at_cell = kwds.pop('node_at_cell', None)
        nodes_at_face = kwds.pop('nodes_at_face', None)

        self._node_at_cell = as_id_array(node_at_cell)
        self._nodes_at_face = as_id_array(nodes_at_face)

        super(DualGraph, self).__init__(*args, **kwds)

        _sort_dual_graph(self)

    @property
    def dual(self):
        return self._dual

    @property
    def node_at_cell(self):
        return self._node_at_cell

    @property
    def cell_at_node(self):
        try:
            return self._cell_at_node
        except AttributeError:
            self._cell_at_node = reverse_one_to_one(
                self.node_at_cell, minlength=self.number_of_nodes)
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
            pair.sort()
            link_at_nodes[tuple(pair)] = link

        link_at_face = np.full((self.number_of_faces, ), -1, dtype=int)
        for face, pair in enumerate(self._nodes_at_face):
            pair.sort()
            link_at_face[face] = link_at_nodes[tuple(pair)]
        self._link_at_face = link_at_face
        return self._link_at_face

    @property
    def face_at_link(self):
        try:
            return self._face_at_link
        except AttributeError:
            self._face_at_link = reverse_one_to_one(
                self.link_at_face, minlength=self.number_of_links)
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

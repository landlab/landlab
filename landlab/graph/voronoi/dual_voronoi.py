from scipy.spatial import Voronoi
import numpy as np

from ...utils.jaggedarray import JaggedArray
from ..dual import DualGraph
from ..graph import Graph
from .voronoi import DelaunayGraph
from .voronoi_helpers import VoronoiConverter
from .voronoi_to_graph import VoronoiDelaunayToGraph


def ugrid_from_voronoi_dual(node_y_and_x, min_cell_size=3,
                            max_node_spacing=None,
                            boundary_nodes=None):
    # from .voronoi import ugrid_from_voronoi
    from ..ugrid import (ugrid_from_unstructured,
                         update_node_at_cell, update_nodes_at_face)

    voronoi = Voronoi(list(zip(node_y_and_x[1], node_y_and_x[0])))

    converter = VoronoiConverter(voronoi, min_patch_size=min_cell_size,
                                 boundary_nodes=boundary_nodes)

    corners = converter.nodes_at_link[:,::-1]
    corners = converter.get_nodes()
    corners = (corners[:, 1], corners[:, 0])
    faces = converter.nodes_at_link
    cells = converter.links_at_patch
    # if len(cells) > 0:
    #     cells = [cell for cell in JaggedArray(cells) if -1 not in cell]

    node_at_cell = converter.get_corner_at_patch()
    nodes_at_face = converter.get_corners_at_link()

    dual = ugrid_from_unstructured(corners, faces, cells)

    return dual, node_at_cell, nodes_at_face


def create_dual_graph(node_y_and_x, min_cell_size=3, max_node_spacing=None):
    voronoi = Voronoi(list(zip(node_y_and_x[1], node_y_and_x[0])))

    converter = VoronoiConverter(voronoi, min_patch_size=min_cell_size)

    corners = converter.get_nodes()
    corners = (corners[:, 1], corners[:, 0])
    faces = converter.get_nodes_at_link()
    cells = converter.get_links_at_patch()
    cells = [cell for cell in JaggedArray(cells)]

    node_at_cell = converter.get_corner_at_patch()
    nodes_at_face = converter.get_corners_at_link()

    graph = Graph(corners, links=faces, patches=cells, sort=False)

    return graph, node_at_cell, nodes_at_face


class DualVoronoiGraph(DualGraph, DelaunayGraph):

    # def __init__(self, node_y_and_x, max_node_spacing=None, min_cell_size=3):
    def __init__(self, node_y_and_x, max_node_spacing=None, sort=False, perimeter_links=None):
        """Create a voronoi grid.

        Parameters
        ----------
        nodes : tuple of array_like
            Coordinates of every node. First *y*, then *x*.

        Examples
        --------
        >>> from landlab.graph import DualVoronoiGraph
        >>> node_x = [0, 1, 2, 3,
        ...           0.2, 1.2, 2.2, 3.2,
        ...           0.4, 1.4, 2.4, 3.4]
        >>> node_y = [0, 0, 0, 0,
        ...           1, 1, 1, 1,
        ...           2, 2, 2, 2]
        >>> graph = DualVoronoiGraph((node_y, node_x), sort=True)
        >>> graph.x_of_corner
        array([ 0.5,  1.5,  2.5,  0.7,  1.7,  2.7,  0.7,  1.7,  2.7,  0.9,  1.9,
                2.9])
        >>> graph.y_of_corner # doctest: +NORMALIZE_WHITESPACE
        array([ 0.42,  0.42,  0.42,  0.58,  0.58,  0.58,  1.42,  1.42,  1.42,
                1.58,  1.58,  1.58])
        >>> graph.corners_at_face # doctest: +NORMALIZE_WHITESPACE
        array([[ 0,  3], [ 3,  1], [ 1,  4], [ 4,  2], [ 2,  5],
               [ 3,  6], [ 4,  7], [ 5,  8],
               [ 6,  9], [ 9,  7], [ 7, 10], [10,  8], [ 8, 11]])
        >>> graph.faces_at_corner # doctest: +NORMALIZE_WHITESPACE
        array([[ 0, -1, -1], [ 2,  1, -1], [ 4,  3, -1],
               [ 5,  0,  1], [ 6,  2,  3], [ 7,  4, -1],
               [ 8,  5, -1], [10,  9,  6], [12, 11,  7],
               [ 8,  9, -1], [10, 11, -1], [12, -1, -1]])
        >>> graph.node_at_cell
        array([5, 6])
        """
        mesh = VoronoiDelaunayToGraph(
            np.vstack((node_y_and_x[1], node_y_and_x[0])).T,
            perimeter_links=perimeter_links,
            # perimeter_links=kwds.get("perimeter_links", None),
        )

        Graph.__init__(
            self,
            node_y_and_x,
            links=mesh.nodes_at_link,
            patches=mesh.links_at_patch,
            sort=False,
        )
        # DelaunayGraph.__init__(
        #     self,
        #     node_y_and_x,
        #     max_node_spacing=max_node_spacing,
        #     sort=False,
        # )

        dual_graph = Graph(
            (mesh.y_of_corner, mesh.x_of_corner),
            links=mesh.corners_at_face,
            patches=mesh.faces_at_cell,
            # node_y_and_x,
            # links=mesh.nodes_at_link,
            # patches=mesh.links_at_patch,
            sort=False,
        )

        # dual_graph, node_at_cell, nodes_at_face = create_dual_graph(
        #     node_y_and_x, min_cell_size=min_cell_size, max_node_spacing=max_node_spacing
        # )

        self.merge(dual_graph, node_at_cell=mesh.node_at_cell, nodes_at_face=mesh.nodes_at_face)

        if sort:
            self.sort()



def __REMOVE_THIS():
        return

        max_node_spacing = kwds.pop('max_node_spacing', None)
        VoronoiGraph.__init__(self, node_y_and_x,
                              max_node_spacing=max_node_spacing, sort=False)

        is_boundary_link = self.patches_at_link[:, 1] == -1

        (dual,
         node_at_cell,
         nodes_at_face) = ugrid_from_voronoi_dual(
             node_y_and_x,
             min_cell_size=min_cell_size,
             max_node_spacing=max_node_spacing,
             boundary_nodes=self.nodes_at_link[is_boundary_link])

        if len(node_at_cell) == 0:
            node_at_cell = None
        self._dual = Graph(dual, sort=False)
        DualGraph.__init__(self, node_at_cell=node_at_cell,
                           nodes_at_face=nodes_at_face, sort=True)
                           # nodes_at_face=nodes_at_face, sort=False)

        return

        max_node_spacing = kwds.pop('max_node_spacing', None)
        (dual,
         node_at_cell,
         nodes_at_face) = ugrid_from_voronoi_dual(node_y_and_x,
                                                  min_cell_size=min_cell_size,
                                                  max_node_spacing=max_node_spacing)
        if len(node_at_cell) == 0:
            node_at_cell = None

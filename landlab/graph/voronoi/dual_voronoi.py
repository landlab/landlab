from scipy.spatial import Voronoi

from ...utils.jaggedarray import JaggedArray
from ..dual import DualGraph
from ..graph import Graph
from .voronoi import VoronoiGraph
from .voronoi_helpers import VoronoiConverter


def ugrid_from_voronoi_dual(node_y_and_x, min_cell_size=3, max_node_spacing=None):
    from ..ugrid import ugrid_from_unstructured

    voronoi = Voronoi(list(zip(node_y_and_x[1], node_y_and_x[0])))

    converter = VoronoiConverter(voronoi, min_patch_size=min_cell_size)

    corners = converter.get_nodes()
    corners = (corners[:, 1], corners[:, 0])
    faces = converter.get_nodes_at_link()
    cells = converter.get_links_at_patch()
    cells = [cell for cell in JaggedArray(cells)]

    node_at_cell = converter.get_corner_at_patch()
    nodes_at_face = converter.get_corners_at_link()

    dual = ugrid_from_unstructured(corners, faces, cells)

    return dual, node_at_cell, nodes_at_face


class DualVoronoiGraph(VoronoiGraph, DualGraph):

    """Dual graph of a voronoi grid."""

    def __init__(self, node_y_and_x, min_cell_size=3, **kwds):
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
        >>> graph = DualVoronoiGraph((node_y, node_x))
        >>> graph.x_of_corner
        array([ 1.5,  2.5,  0.7,  1.7,  2.7,  0.7,  1.7,  2.7,  0.9,  1.9])
        >>> graph.y_of_corner # doctest: +NORMALIZE_WHITESPACE
        array([ 0.42,  0.42,  0.58,  0.58,  0.58,  1.42,  1.42,  1.42,  1.58,
                1.58])
        >>> graph.corners_at_face # doctest: +NORMALIZE_WHITESPACE
        array([[2, 0], [0, 3], [3, 1], [1, 4],
               [2, 5], [3, 6], [4, 7],
               [5, 8], [8, 6], [6, 9], [9, 7]])
        >>> graph.faces_at_corner # doctest: +NORMALIZE_WHITESPACE
        array([[ 1,  0, -1], [ 3,  2, -1],
               [ 4,  0, -1], [ 5,  1,  2], [ 6,  3, -1],
               [ 7,  4, -1], [ 9,  8,  5], [10,  6, -1],
               [ 7,  8, -1], [ 9, 10, -1]])
        >>> graph.node_at_cell
        array([5, 6])
        """
        max_node_spacing = kwds.pop("max_node_spacing", None)
        (dual, node_at_cell, nodes_at_face) = ugrid_from_voronoi_dual(
            node_y_and_x, min_cell_size=min_cell_size, max_node_spacing=max_node_spacing
        )

        self._dual = Graph(dual, sort=False)
        VoronoiGraph.__init__(
            self, node_y_and_x, max_node_spacing=max_node_spacing, sort=False
        )
        DualGraph.__init__(
            self, node_at_cell=node_at_cell, nodes_at_face=nodes_at_face, sort=True
        )

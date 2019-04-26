import numpy as np

from ..dual import DualGraph
from .structured_quad import StructuredQuadGraph


class DualStructuredQuadGraphExtras(object):
    @property
    def corners_at_right_edge(self):
        return self.dual.nodes_at_right_edge

    @property
    def corners_at_top_edge(self):
        return self.dual.nodes_at_top_edge

    @property
    def corners_at_left_edge(self):
        return self.dual.nodes_at_left_edge

    @property
    def corners_at_bottom_edge(self):
        return self.dual.nodes_at_bottom_edge

    @property
    def perimeter_corners(self):
        return self.dual.perimeter_nodes

    @property
    def horizontal_faces(self):
        return self.dual.horizontal_links

    @property
    def vertical_faces(self):
        return self.dual.vertical_links


class DualStructuredQuadGraph(
    DualStructuredQuadGraphExtras, StructuredQuadGraph, DualGraph
):

    """Dual graph of a structured grid of quadrilaterals.

    Examples
    --------
    >>> from landlab.graph import DualStructuredQuadGraph
    >>> node_y = [-1, -2, -3,
    ...            0,  0,  0,
    ...            1,  2,  3]
    >>> node_x = [ 0,  1,  2,
    ...            0,  2,  3,
    ...            0,  1,  2]
    >>> graph = DualStructuredQuadGraph((node_y, node_x), shape=(3, 3))
    >>> graph.number_of_corners == 4
    True
    >>> graph.y_of_corner
    array([-1.25, -0.75,  0.75,  1.25])
    >>> graph.x_of_corner
    array([ 2.  ,  0.75,  0.75,  2.  ])
    >>> graph.node_at_cell
    array([4])
    """

    def __init__(self, node_y_and_x, shape=None):
        dual_y, dual_x = get_corners(node_y_and_x, shape)
        dual_shape = dual_y.shape
        node_at_cell = get_node_at_cell(shape)
        nodes_at_face = get_nodes_at_face(shape)

        self._dual = StructuredQuadGraph((dual_y, dual_x), shape=dual_shape)
        StructuredQuadGraph.__init__(self, node_y_and_x, shape=shape)
        DualGraph.__init__(self, node_at_cell=node_at_cell, nodes_at_face=nodes_at_face)


class DualRectilinearGraph(DualStructuredQuadGraph):

    """Create a dual graph for a rectilinear grid.

    Examples
    --------
    >>> from landlab.graph import DualRectilinearGraph
    >>> graph = DualRectilinearGraph(([0, 1, 3], [0, 5, 15, 30]))
    >>> graph.x_of_corner # doctest: +NORMALIZE_WHITESPACE
    array([ 2.5, 10. , 22.5,
            2.5, 10. , 22.5])
    >>> graph.y_of_corner # doctest: +NORMALIZE_WHITESPACE
    array([ 0.5, 0.5, 0.5,
            2. , 2. , 2. ])
    >>> graph.number_of_cells == 2
    True
    >>> graph.faces_at_cell
    array([[3, 5, 2, 0],
           [4, 6, 3, 1]])
    """

    def __init__(self, node_y_and_x):
        shape = (len(node_y_and_x[0]), len(node_y_and_x[1]))
        node_y_and_x = np.meshgrid(*node_y_and_x, indexing="ij")

        super(DualRectilinearGraph, self).__init__(node_y_and_x, shape)


class DualUniformRectilinearGraph(DualRectilinearGraph):

    """Create a dual graph for a uniform rectilinear grid.

    Examples
    --------
    >>> from landlab.graph import DualUniformRectilinearGraph
    >>> graph = DualUniformRectilinearGraph((4, 3))
    >>> graph.x_of_corner # doctest: +NORMALIZE_WHITESPACE
    array([ 0.5, 1.5,
            0.5, 1.5,
            0.5, 1.5])
    >>> graph.y_of_corner # doctest: +NORMALIZE_WHITESPACE
    array([ 0.5, 0.5,
            1.5, 1.5,
            2.5, 2.5])
    >>> graph.number_of_cells == 2
    True
    >>> graph.faces_at_cell
    array([[2, 3, 1, 0],
           [5, 6, 4, 3]])
    """

    def __init__(self, shape, spacing=(1.0, 1.0), origin=(0.0, 0.0)):
        spacing = np.broadcast_to(spacing, 2)
        origin = np.broadcast_to(origin, 2)

        node_y_and_x = (
            np.arange(shape[0]) * spacing[0] + origin[0],
            np.arange(shape[1]) * spacing[1] + origin[1],
        )

        super(DualUniformRectilinearGraph, self).__init__(node_y_and_x)


def get_node_at_cell(shape):
    """Set up an array that gives the node at each cell.

    Examples
    --------
    >>> from landlab.graph.structured_quad.dual_structured_quad import get_node_at_cell
    >>> get_node_at_cell((5, 6)) # doctest: +NORMALIZE_WHITESPACE
    array([ 7,  8,  9, 10,
           13, 14, 15, 16,
           19, 20, 21, 22])
    """
    from .ext.at_cell import fill_node_at_cell

    node_at_cell = np.empty((shape[0] - 2) * (shape[1] - 2), dtype=int)

    fill_node_at_cell(shape, node_at_cell)

    return node_at_cell


def get_nodes_at_face(shape):
    """Set up an array that gives the nodes on either side of each face.

    Examples
    --------
    >>> from landlab.graph.structured_quad.dual_structured_quad import get_nodes_at_face
    >>> get_nodes_at_face((3, 4)) # doctest: +NORMALIZE_WHITESPACE
    array([[ 1,  5], [ 2,  6],
           [ 4,  5], [ 5,  6], [ 6,  7],
           [ 5,  9], [ 6, 10]])
    """
    from .ext.at_face import fill_nodes_at_face

    n_faces = (shape[1] - 2) * (shape[0] - 1) + (shape[0] - 2) * (shape[1] - 1)
    nodes_at_face = np.empty((n_faces, 2), dtype=int)
    fill_nodes_at_face(shape, nodes_at_face)

    return nodes_at_face


def get_corners(node_y_and_x, shape):
    y_of_node, x_of_node = (
        np.asarray(node_y_and_x[0], dtype=float),
        np.asarray(node_y_and_x[1], dtype=float),
    )
    y_of_node.shape = x_of_node.shape = shape

    x_of_corner = (
        x_of_node[:-1, :-1]
        + x_of_node[:-1, 1:]
        + x_of_node[1:, :-1]
        + x_of_node[1:, 1:]
    ) * 0.25
    y_of_corner = (
        y_of_node[:-1, :-1]
        + y_of_node[:-1, 1:]
        + y_of_node[1:, :-1]
        + y_of_node[1:, 1:]
    ) * 0.25

    return y_of_corner, x_of_corner


def reshape_nodes(y_and_x, shape=None):
    y_of_node, x_of_node = (
        np.asarray(y_and_x[0], dtype=float),
        np.asarray(y_and_x[1], dtype=float),
    )

    if y_of_node.size != x_of_node.size:
        raise ValueError("size mismatch for size of x and y")

    if shape is None:
        if y_of_node.shape != x_of_node.shape:
            raise ValueError("shape mismatch for size of x and y")
    else:
        y_of_node.shape = shape
        x_of_node.shape = shape
    return y_of_node, x_of_node

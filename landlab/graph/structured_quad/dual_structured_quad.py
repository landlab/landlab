import numpy as np

from ..dual import DualGraph
from .structured_quad import (StructuredQuadGraph, RectilinearGraph,
                              UniformRectilinearGraph, )


class DualStructuredQuadGraph(DualGraph, StructuredQuadGraph):

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
    >>> graph.number_of_corners
    4
    >>> graph.y_of_corner
    array([-1.25, -0.75,  0.75,  1.25])
    >>> graph.x_of_corner
    array([ 2.  ,  0.75,  0.75,  2.  ])
    >>> graph.node_at_cell
    array([4])
    """

    def __init__(self, nodes, shape=None):
        y_of_node, x_of_node = reshape_nodes(nodes, shape=shape)
        shape = y_of_node.shape

        dual_y, dual_x = get_corners((y_of_node, x_of_node), shape)

        dual_shape = dual_y.shape

        self._dual = StructuredQuadGraph((dual_y, dual_x), shape=dual_shape)

        node_at_cell = get_node_at_cell(shape)
        nodes_at_face = get_nodes_at_face(shape)

        super(DualStructuredQuadGraph, self).__init__(
            (y_of_node, x_of_node), shape=shape, node_at_cell=node_at_cell,
            nodes_at_face=nodes_at_face)

    @property
    def perimeter_corners(self):
        return self.dual.perimeter_nodes


class DualRectilinearGraph(DualGraph, RectilinearGraph):

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
    >>> graph.number_of_cells
    2
    >>> graph.faces_at_cell
    array([[3, 5, 2, 0],
           [4, 6, 3, 1]])
    """

    def __init__(self, nodes):
        dual_nodes = (np.asarray(nodes[0], dtype=float),
                      np.asarray(nodes[1], dtype=float))
        dual_nodes = [x[:-1] + np.diff(x) * .5 for x in dual_nodes]

        self._dual = RectilinearGraph(dual_nodes)
        # self._node_at_cell = get_node_at_cell(self.shape)

        shape = (len(nodes[0]), len(nodes[1]))

        node_at_cell = get_node_at_cell(shape)
        nodes_at_face = get_nodes_at_face(shape)

        super(DualRectilinearGraph, self).__init__(nodes,
                                                   node_at_cell=node_at_cell,
                                                   nodes_at_face=nodes_at_face)


class DualUniformRectilinearGraph(DualGraph, UniformRectilinearGraph):

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
    >>> graph.number_of_cells
    2
    >>> graph.faces_at_cell
    array([[2, 3, 1, 0],
           [5, 6, 4, 3]])
    """

    def __init__(self, shape, spacing=(1., 1.), origin=(0., 0.)):
        spacing = np.broadcast_to(spacing, 2)
        origin = np.broadcast_to(origin, 2)

        dual_shape = [dim - 1 for dim in shape]
        dual_origin = [x + dx * .5 for x, dx in zip(origin, spacing)]

        self._dual = UniformRectilinearGraph(dual_shape,
                                             spacing=spacing,
                                             origin=dual_origin)

        node_at_cell = get_node_at_cell(shape)
        nodes_at_face = get_nodes_at_face(shape)

        super(DualUniformRectilinearGraph, self).__init__(
            shape, spacing=spacing, origin=origin, node_at_cell=node_at_cell,
            nodes_at_face=nodes_at_face)


    @property
    def dual(self):
        return self._dual


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


def get_corners(nodes, shape):
    y_of_node, x_of_node = (np.array(nodes[0], dtype=float),
                            np.array(nodes[1], dtype=float))
    if shape is not None:
        y_of_node.shape = shape
        x_of_node.shape = shape

    x_of_corner = (x_of_node[:-1, :-1] + x_of_node[:-1, 1:] +
                   x_of_node[1:, :-1] + x_of_node[1:, 1:]) * .25
    y_of_corner = (y_of_node[:-1, :-1] + y_of_node[:-1, 1:] +
                   y_of_node[1:, :-1] + y_of_node[1:, 1:]) * .25

    return y_of_corner, x_of_corner


def reshape_nodes(nodes, shape=None):
    y_of_node, x_of_node = (np.asarray(nodes[0], dtype=float),
                            np.asarray(nodes[1], dtype=float))

    if y_of_node.size != x_of_node.size:
        raise ValueError('size mismatch for size of x and y')

    if shape is None:
        if y_of_node.shape != x_of_node.shape:
            raise ValueError('shape mismatch for size of x and y')
    else:
        y_of_node.shape = shape
        x_of_node.shape = shape
    return y_of_node, x_of_node

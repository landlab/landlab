import numpy as np

from ..dual import DualGraph
from .structured_quad import RectilinearGraph
from .structured_quad import StructuredQuadGraph
from .structured_quad import UniformRectilinearGraph


class DualStructuredQuadGraph(DualGraph, StructuredQuadGraph):
    """Dual graph of a structured grid of quadrilaterals.

    Examples
    --------
    >>> from landlab.graph import DualStructuredQuadGraph
    >>> node_y = [-1, -2, -3, 0, 0, 0, 1, 2, 3]
    >>> node_x = [0, 1, 2, 0, 2, 3, 0, 1, 2]
    >>> graph = DualStructuredQuadGraph((node_y, node_x), shape=(3, 3), sort=True)
    >>> graph.number_of_corners == 4
    True
    >>> graph.y_of_corner
    array([-1.25, -0.75,  0.75,  1.25])
    >>> graph.x_of_corner
    array([2.  , 0.75, 0.75, 2.  ])
    >>> graph.node_at_cell
    array([4])
    """

    def __init__(self, node_y_and_x, shape=None, sort=True):
        StructuredQuadGraph.__init__(self, node_y_and_x, shape=shape)

        dual_graph = StructuredQuadGraph(
            DualStructuredQuadGraph.get_corners(node_y_and_x, self.shape)
        )

        self.merge(
            dual_graph,
            node_at_cell=DualStructuredQuadGraph.get_node_at_cell(self.shape),
            nodes_at_face=DualStructuredQuadGraph.get_nodes_at_face(self.shape),
        )

        if sort:
            self.sort()

    @staticmethod
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

    @staticmethod
    def get_node_at_cell(shape):
        """Set up an array that gives the node at each cell.

        Examples
        --------
        >>> from landlab.graph.structured_quad import DualStructuredQuadGraph
        >>> DualStructuredQuadGraph.get_node_at_cell((5, 6))
        array([ 7,  8,  9, 10,
               13, 14, 15, 16,
               19, 20, 21, 22])
        """
        from .ext.at_cell import fill_node_at_cell

        node_at_cell = np.empty((shape[0] - 2) * (shape[1] - 2), dtype=int)

        fill_node_at_cell(shape, node_at_cell)

        return node_at_cell

    @staticmethod
    def get_nodes_at_face(shape):
        """Set up an array that gives the nodes on either side of each face.

        Examples
        --------
        >>> from landlab.graph.structured_quad import DualStructuredQuadGraph
        >>> DualStructuredQuadGraph.get_nodes_at_face((3, 4))
        array([[ 1,  5], [ 2,  6],
               [ 4,  5], [ 5,  6], [ 6,  7],
               [ 5,  9], [ 6, 10]])
        """
        from .ext.at_face import fill_nodes_at_face

        n_faces = (shape[1] - 2) * (shape[0] - 1) + (shape[0] - 2) * (shape[1] - 1)
        nodes_at_face = np.empty((n_faces, 2), dtype=int)
        fill_nodes_at_face(shape, nodes_at_face)

        return nodes_at_face


class DualRectilinearGraph(DualGraph, RectilinearGraph):
    """Create a dual graph for a rectilinear grid.

    Examples
    --------
    >>> from landlab.graph import DualRectilinearGraph
    >>> graph = DualRectilinearGraph(([0, 1, 3], [0, 5, 15, 30]))
    >>> graph.x_of_corner.reshape((2, 3))
    array([[  2.5,  10. ,  22.5],
           [  2.5,  10. ,  22.5]])
    >>> graph.y_of_corner.reshape((2, 3))
    array([[0.5, 0.5, 0.5],
           [2. , 2. , 2. ]])
    >>> graph.number_of_cells == 2
    True
    >>> graph.faces_at_cell
    array([[3, 5, 2, 0],
           [4, 6, 3, 1]])
    """

    def __init__(self, node_y_and_x):
        RectilinearGraph.__init__(self, node_y_and_x)

        dual_graph = RectilinearGraph(DualRectilinearGraph.get_corners(node_y_and_x))

        self.merge(
            dual_graph,
            node_at_cell=DualStructuredQuadGraph.get_node_at_cell(self.shape),
            nodes_at_face=DualStructuredQuadGraph.get_nodes_at_face(self.shape),
        )

    @staticmethod
    def get_corners(node_y_and_x):
        y_of_node, x_of_node = (
            np.asarray(node_y_and_x[0]),
            np.asarray(node_y_and_x[1]),
        )

        return (
            (y_of_node[1:] + y_of_node[:-1]) * 0.5,
            (x_of_node[1:] + x_of_node[:-1]) * 0.5,
        )


class DualUniformRectilinearGraph(DualGraph, UniformRectilinearGraph):
    """Create a dual graph for a uniform rectilinear grid.

    Examples
    --------
    >>> from landlab.graph import DualUniformRectilinearGraph
    >>> graph = DualUniformRectilinearGraph((4, 3))
    >>> graph.x_of_corner.reshape((3, 2))
    array([[0.5, 1.5],
           [0.5, 1.5],
           [0.5, 1.5]])
    >>> graph.y_of_corner.reshape((3, 2))
    array([[0.5, 0.5],
           [1.5, 1.5],
           [2.5, 2.5]])
    >>> graph.number_of_cells == 2
    True
    >>> graph.faces_at_cell
    array([[2, 3, 1, 0],
           [5, 6, 4, 3]])
    """

    def __init__(self, shape, spacing=1.0, origin=(0.0, 0.0)):
        spacing = np.broadcast_to(spacing, 2)
        origin = np.broadcast_to(origin, 2)

        UniformRectilinearGraph.__init__(self, shape, spacing=spacing, origin=origin)

        dual_graph = UniformRectilinearGraph(
            (shape[0] - 1, shape[1] - 1), spacing=spacing, origin=origin + spacing * 0.5
        )

        self.merge(
            dual_graph,
            node_at_cell=DualStructuredQuadGraph.get_node_at_cell(self.shape),
            nodes_at_face=DualStructuredQuadGraph.get_nodes_at_face(self.shape),
        )

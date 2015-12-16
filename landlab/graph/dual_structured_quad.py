import numpy as np

from .dual import DualGraphMixIn
from .structured_quad import (StructuredQuadGraph, RectilinearGraph,
                              UniformRectilinearGraph, )


class DualStructuredQuadGraph(StructuredQuadGraph, DualGraphMixIn):

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
    array([-0.75, -1.25,  0.75,  1.25])
    >>> graph.x_of_corner
    array([ 0.75,  2.  ,  0.75,  2.  ])
    """

    def __init__(self, nodes, shape=None):
        super(DualStructuredQuadGraph, self).__init__(nodes, shape=shape)

        dual_y = np.mean(self.y_of_node[self.nodes_at_patch], axis=1)
        dual_x = np.mean(self.x_of_node[self.nodes_at_patch], axis=1)

        dual_shape = [dim - 1 for dim in self.shape]

        self._dual = StructuredQuadGraph((dual_y, dual_x), shape=dual_shape)


class DualRectilinearGraph(RectilinearGraph, DualGraphMixIn):

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
        super(DualRectilinearGraph, self).__init__(nodes)

        dual_nodes = (np.asarray(nodes[0], dtype=float),
                      np.asarray(nodes[1], dtype=float))
        dual_nodes = [x[:-1] + np.diff(x) * .5 for x in dual_nodes]

        self._dual = RectilinearGraph(dual_nodes)


class DualUniformRectilinearGraph(UniformRectilinearGraph, DualGraphMixIn):

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
        super(DualUniformRectilinearGraph, self).__init__(shape,
                                                          spacing=spacing,
                                                          origin=origin)

        dual_shape = [dim - 1 for dim in shape]
        dual_origin = [x + dx * .5 for x, dx in zip(origin, spacing)]

        self._dual = UniformRectilinearGraph(dual_shape,
                                             spacing=spacing,
                                             origin=dual_origin)

    @property
    def dual(self):
        return self._dual

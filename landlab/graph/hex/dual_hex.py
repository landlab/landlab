import numpy as np

from ..voronoi import DualVoronoiGraph
from .hex import create_xy_of_node


class DualHexGraph(DualVoronoiGraph):

    """Graph of a structured grid of triangles.

    Examples
    --------
    >>> from landlab.graph import HexGraph
    >>> graph = StructuredQuadGraph((3, 2))
    >>> graph.number_of_nodes
    7
    >>> graph.y_of_node # doctest: +NORMALIZE_WHITESPACE
    array([ 0.,  0.,  1.,  1.,  1.,  2.,  2.])
    >>> graph.x_of_node # doctest: +NORMALIZE_WHITESPACE
    array([ 0.5,  1.5,  0. ,  1. ,  2. ,  0.5,  1.5])
    """

    def __init__(self, shape, spacing=1., origin=(0., 0.),
                 orientation='horizontal'):
        """Create a structured grid of triangles.

        Parameters
        ----------
        shape : tuple of int
            Number of rows and columns of the hex grid. The first value
            is the number of nodes in the first column and the second the
            number of nodes in the first column.
        spacing : float, optional
            Length of links.
        origin : tuple of float, optional
            Coordinates of lower-left corner of the grid.
        """
        try:
            spacing = float(spacing)
        except TypeError:
            raise TypeError('spacing must be a float')

        x_of_node, y_of_node = create_xy_of_node(shape, spacing=spacing,
                                                 origin=origin,
                                                 orientation=orientation)

        super(DualHexGraph, self).__init__(
            (y_of_node, x_of_node), xy_sort=True, rot_sort=True,
            min_cell_size=6, max_node_spacing=shape[1] + 1)

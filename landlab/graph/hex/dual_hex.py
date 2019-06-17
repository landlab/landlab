from ..voronoi import DualVoronoiGraph
from .hex import setup_xy_of_node


class DualHexGraph(DualVoronoiGraph):

    """Graph of a structured grid of triangles.

    Examples
    --------
    >>> import numpy as np
    >>> from landlab.graph import HexGraph

    >>> graph = DualHexGraph((3, 2), node_layout='hex')
    >>> graph.number_of_nodes
    7
    >>> graph.number_of_corners
    6

    >>> np.round(graph.y_of_node * 2. / np.sqrt(3))
    ...     # doctest: +NORMALIZE_WHITESPACE
    array([ 0.,  0.,  1.,  1.,  1.,  2.,  2.])
    >>> graph.x_of_node # doctest: +NORMALIZE_WHITESPACE
    array([ 0. ,  1. , -0.5,  0.5,  1.5,  0. ,  1. ])
    """

    def __init__(
        self,
        shape,
        spacing=1.0,
        origin=(0.0, 0.0),
        orientation="horizontal",
        node_layout="rect",
    ):
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
        orientation: {'horizontal', 'vertical'}
            Specify if triangles should be laid out in rows or columns.
        node_layout: {'rect', 'hex', 'rect1'}
            Specify the overall layout of the nodes. Use *rect* for
            the layout to approximate a rectangle and *hex* for
            a hexagon.
        """
        try:
            spacing = float(spacing)
        except TypeError:
            raise TypeError("spacing must be a float")

        x_of_node, y_of_node = setup_xy_of_node(
            shape,
            spacing=spacing,
            origin=origin,
            orientation=orientation,
            node_layout=node_layout,
        )

        if node_layout == "hex":
            max_node_spacing = None
        else:
            max_node_spacing = shape[1] + 1

        DualVoronoiGraph.__init__(
            self,
            (y_of_node, x_of_node),
            xy_sort=True,
            rot_sort=True,
            min_cell_size=6,
            max_node_spacing=max_node_spacing,
        )
